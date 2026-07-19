#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##################### DIFFERENTIAL GENE EXPRESSION ANALYSIS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Description:          Bulletproof DGE Analysis with nested tryCatch and 
#                       cell count pre-validation to prevent ANY Seurat crash.

library(Seurat)
library(dplyr)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
rds_dir <- file.path(base_dir, "rds_OUT")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  input_filename <- "IntgratedRO_LogNorm_1stLvl.rds"
} else {
  input_filename <- args[1]
}

input_path <- file.path(rds_dir, input_filename)

if (!file.exists(input_path)) {
  stop("인풋 파일이 존재하지 않습니다: ", input_path)
}

message("Loading Seurat object...")
data <- readRDS(input_path)

DefaultAssay(data) <- "RNA"
group <- "disease"
Idents(data) <- "seurat_clusters"
top <- 500

message("Normalizing Data...")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

marker <- list()

# 모든 클러스터에 대해 루프를 돕니다.
for (cluster in levels(data$seurat_clusters)) {
  message("\nProcessing cluster: ", cluster)
  
  # 1. 사전 검증: 클러스터의 총 세포 수 확인
  cell_count <- length(WhichCells(data, idents = cluster))
  if (cell_count < 3) {
    message("  -> [건너뜀] 클러스터 ", cluster, " 의 세포 수가 너무 적습니다 (", cell_count, "개).")
    next
  }
  
  # 2. 무적의 중첩 tryCatch 블록
  x <- tryCatch({
    message("  -> FindConservedMarkers 시도 중...")
    res <- FindConservedMarkers(data, assay = "RNA", layer = "data",
                                ident.1 = cluster, grouping.var = group,
                                only.pos = TRUE)
    
    # 결과가 존재하고 최소 p-value 컬럼이 있을 때만 필터링
    if (!is.null(res) && nrow(res) > 0 && "minimump_p_val" %in% colnames(res)) {
      res <- res %>% slice_min(minimump_p_val, n = top)
    }
    res
    
  }, error = function(e) {
    message("  -> [우회 동작] 한쪽 그룹 세포 부족 등으로 인해 FindMarkers로 자동 전환합니다.")
    
    tryCatch({
      res2 <- FindMarkers(data, layer = "data", ident.1 = cluster, only.pos = TRUE)
      
      # 결과가 존재하고 조정된 p-value 컬럼이 있을 때만 필터링
      if (!is.null(res2) && nrow(res2) > 0 && "p_val_adj" %in% colnames(res2)) {
        res2 <- res2 %>% slice_min(p_val_adj, n = top)
      }
      res2
      
    }, error = function(e2) {
      message("  -> [완전 실패] 해당 클러스터는 마커 탐색 조건에 부합하지 않아 건너뜁니다.")
      return(NULL)
    })
  })
  
  # 결과물을 리스트에 저장
  if (!is.null(x) && nrow(x) > 0) {
    x$cluster <- cluster
    marker[[length(marker) + 1]] <- x
  }
}

if (!dir.exists(rds_dir)) {
  dir.create(rds_dir, recursive = TRUE)
}

output_file <- file.path(rds_dir, "DEGs_IntegratedRO_LogNorm_1stLvl.rds")
saveRDS(marker, file = output_file, compress = TRUE)

cat("\n[완료] DGE 분석이 모든 악조건을 뚫고 끝까지 완료되었습니다!\n마커 파일 저장 경로: ", output_file, "\n")