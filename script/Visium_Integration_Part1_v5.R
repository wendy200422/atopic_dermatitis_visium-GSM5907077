#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
################ NORMALIZATION AND INTEGRATION OF VISIUM DATASET ###############
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Description:          Complete integration pipeline - QC, Filtering, Normalization.

library(Seurat)
library(dplyr)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
setwd(base_dir)
rds_dir <- file.path(base_dir, "rds_OUT")

if (!dir.exists(rds_dir)) {
  dir.create(rds_dir, recursive = TRUE)
}

# [샘플 목록 확보] GSM5907096 폴더 하위의 모든 폴더를 탐색
# 'rds_OUT', 'visualization' 등 분석 폴더를 제외하고 오직 샘플 폴더만 필터링합니다.
all_folders <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)
visium_samples <- all_folders[!grepl("rds_OUT|visualization", all_folders)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############## II. PREPROCESSING: QC AND NORMALIZATION ########################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

stdata_visium <- list()

for (i in 1:length(visium_samples)) {
  id <- basename(visium_samples[i])
  message("\n[", i, "/", length(visium_samples), "] Processing: ", id)
  
  # Visium 데이터 로드
  x <- Load10X_Spatial(data.dir = visium_samples[i], slice = id)
  x$orig.ident <- id
  
  # QC Metrics: 미토콘드리아 및 헤모글로빈 비율 계산
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x[["percent.hgb"]] <- PercentageFeatureSet(x, features = intersect(c("HBA1", "HBA2", "HBB"), rownames(x)))
  
  # 필터링: MT 20%, Hb 10%, UMI 200 이상, 유전자 100 이상인 스팟만 유지
  x <- subset(x, subset = percent.hgb < 10 & percent.mt < 20 & nCount_Spatial > 200 & nFeature_Spatial > 100)
  
  stdata_visium[[length(stdata_visium) + 1]] <- x
  rm(x) # 메모리 관리
}

# [스팟 개수 검증] 통합 시 에러 방지를 위해 스팟이 100개 미만인 샘플은 자동 제거
stdata_visium <- stdata_visium[sapply(stdata_visium, ncol) >= 100]

# 정규화 및 스케일링
message("\n--- Normalizing and Scaling data ---")
stdata_visium <- lapply(X = stdata_visium, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  return(x)
})

# [저장] 최종 결과물 저장
output_file <- file.path(rds_dir, "SampleListVisium_LogNorm.rds")
saveRDS(stdata_visium, file = output_file, compress = TRUE)

cat("\n[완료] Visium 데이터 전처리가 완벽하게 끝났습니다!\n")
cat("저장 경로: ", output_file, "\n")
cat("성공적으로 처리된 샘플 수: ", length(stdata_visium), "\n")