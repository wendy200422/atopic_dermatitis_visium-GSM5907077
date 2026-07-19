#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
########### IV. DIMENSION REDUCTION AND CLUSTERING ON INTEGRATED DATA ##########
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Description:          Perfect fix for ScaleData multi-layer conflicts, 
#                       S4 casting errors (as.Assay), RunPCA feature mapping errors,
#                       disease metadata parsing, and color palette limits.

# Loads required packages
library(Seurat)
library(SeuratObject) # S4 메서드(as)가 원활히 작동하도록 명시적으로 로드합니다.
library(dplyr)
library(cowplot)
library(rcartocolor)

# [경로 설정] 사용자의 로컬 윈도우 디렉토리 경로
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
rds_dir <- file.path(base_dir, "rds_OUT")

# 2단계 통합 결과 파일 로드 (Cell 유실 없는 rpca 통합 객체)
input_file <- file.path(rds_dir, "IntegratedRO_LogNorm.rds")
if (!file.exists(input_file)) {
  stop("2단계 결과 파일이 존재하지 않습니다. 먼저 Part2를 다시 실행해주세요: ", input_file)
}

message("Loading IntegratedRO_LogNorm.rds...")
scdata_integrated <- readRDS(input_file)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# [Seurat v5 대응 - 클래스 캐스팅 에러 완벽 해결]
# 1. integrated 어세이를 기본으로 설정합니다.
# 2. 다중 레이어 충돌을 방지하기 위해 v5 Assay5 형태를 v3 Assay 형태로 변환합니다.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
DefaultAssay(scdata_integrated) <- "integrated"
scdata_integrated[["integrated"]] <- as(object = scdata_integrated[["integrated"]], Class = "Assay")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# [메타데이터 오류 해결 - 질환(disease) 그룹핑]
# 샘플 이름에 "AD"가 포함되어 있으면 "AD", 그렇지 않으면 "HC"로 자동 분류합니다.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
scdata_integrated$disease <- ifelse(grepl("AD", scdata_integrated$orig.ident, ignore.case = TRUE), "AD", "HC")

# 정상적으로 2개의 그룹(AD, HC)으로 잘 나뉘었는지 콘솔에 출력해서 확인합니다.
message("Disease group counts:")
print(table(scdata_integrated$disease))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# [데이터 스케일링 및 PCA 재연산]
# 단독 Assay 변환 후, 명시적(features)으로 대상을 지정하여 PCA를 재생성합니다.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
message("Scaling data...")
integrated_genes <- rownames(scdata_integrated[["integrated"]])
scdata_integrated <- ScaleData(scdata_integrated, features = integrated_genes)

# [오류 해결] features 인자를 명시하여 PCA 연산 실패(Warning)를 방지합니다.
message("Running PCA (Generating 'pca' DimReduc)...")
scdata_integrated <- RunPCA(scdata_integrated, features = integrated_genes, npcs = 50, verbose = FALSE, reduction.name = "pca")

# Finds neighbors in PCA space (using 35 PCs as per original design)
message("Finding neighbors...")
scdata_integrated <- FindNeighbors(scdata_integrated, reduction = "pca", dims = 1:35)

# Finds clusters using Louvain Algorithm (Resolution = 0.02)
message("Finding clusters...")
scdata_integrated <- FindClusters(scdata_integrated, resolution = 0.02)

# Runs UMAP dimensional reduction
message("Running UMAP...")
scdata_integrated <- RunUMAP(scdata_integrated, reduction = "pca", dims = 1:35, reduction.name = "umap")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# [시각화 오류 해결 - 색상 팔레트 자동 확장]
# 클러스터 수가 12개를 초과하더라도 에러가 발생하지 않도록 그라데이션으로 색상을 생성합니다.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
plot_dir <- file.path(base_dir, "visualization", "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

message("Saving UMAP plots...")
clst.nr <- length(levels(scdata_integrated@meta.data[["seurat_clusters"]]))

# 13개 이상일 때 에러 방지 (n 파라미터 지정)
safe_base_cols <- carto_pal(n = min(clst.nr, 12), name = "Safe") 
clst.cols <- colorRampPalette(safe_base_cols)(clst.nr)

png(file.path(plot_dir, "UMAP_Integrated_Rohjan_LogNorm_1stLvl.png"),
    width = 600, height = 200, units = 'mm', res = 400)

p1 <- DimPlot(scdata_integrated, reduction = "umap", group.by = "disease",
              cols = c("#CC6677", "#88CCEE"), shuffle = TRUE)

p2 <- DimPlot(scdata_integrated, reduction = "umap", label = TRUE,
              repel = TRUE, cols = clst.cols, shuffle = TRUE, label.size = 5)

print(plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(1/2, 1/2)))
dev.off()

# Writes the integrated 1st-level Seurat object to a .rds file
output_file <- file.path(rds_dir, "IntgratedRO_LogNorm_1stLvl.rds")
saveRDS(scdata_integrated, file = output_file, compress = TRUE)

cat("\n[완료] 3단계 차원축소 및 클러스터링이 에러 없이 완벽히 끝났습니다!\n결과 파일 저장 경로: ", output_file, "\n")