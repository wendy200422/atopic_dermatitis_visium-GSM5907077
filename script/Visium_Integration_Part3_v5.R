#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
########### IV. DIMENSION REDUCTION AND CLUSTERING ON INTEGRATED DATA ##########
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(Seurat)
library(dplyr)
library(cowplot)
library(RColorBrewer)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
setwd(base_dir)
rds_dir  <- file.path(base_dir, "rds_OUT")
plot_dir <- file.path(base_dir, "visualization", "plots")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. 통합 데이터 로드
message("1. Loading integrated data...")
visium_integrated <- readRDS(file.path(rds_dir, "IntegratedVisium_LogNorm.rds"))

# 2. 어세이 설정
message("2. Setting Default Assay to 'integrated'...")
DefaultAssay(visium_integrated) <- "integrated"
# (주의: JoinLayers는 객체 구조상 불필요하므로 제거했습니다)

# 3. 질환군(disease) 메타데이터 생성
message("3. Annotating disease groups...")
visium_integrated$disease <- case_when(
  grepl("LS", visium_integrated$orig.ident) ~ "LS",
  grepl("NL", visium_integrated$orig.ident) ~ "NL",
  grepl("HE", visium_integrated$orig.ident) ~ "HC",
  TRUE ~ "Unknown"
)

# 4. 차원 축소 (PCA & UMAP)
message("4. Running Dimensional Reduction (PCA & UMAP)...")

# [핵심 수정] IntegrateData가 만들어둔 유전자 목록을 그대로 가져옵니다.
# FindVariableFeatures를 다시 실행하면 유전자 목록이 깨지므로 절대 실행하지 않습니다.
integ_features <- VariableFeatures(visium_integrated)

# 가져온 유전자 목록으로 명시적 스케일링 및 PCA 수행
visium_integrated <- ScaleData(visium_integrated, features = integ_features, verbose = FALSE)
visium_integrated <- RunPCA(visium_integrated, features = integ_features, npcs = 30, verbose = FALSE)

visium_integrated <- RunUMAP(visium_integrated, dims = 1:30, verbose = FALSE)

# 5. 공간 클러스터링
message("5. Clustering cells/spots...")
visium_integrated <- FindNeighbors(visium_integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
visium_integrated <- FindClusters(visium_integrated, resolution = 0.5, verbose = FALSE)

# 6. UMAP 시각화 저장
message("6. Saving UMAP plots...")
png(file.path(plot_dir, "UMAP_Visium_Integrated.png"), width = 600, height = 300, units = 'mm', res = 400)
p1 <- DimPlot(visium_integrated, reduction = "umap", group.by = "disease", pt.size = 0.5, shuffle = TRUE)
p2 <- DimPlot(visium_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)
print(plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(1/2, 1/2)))
dev.off()

# 7. 공간 클러스터 시각화 (Spatial Plot) 저장
message("7. Saving Spatial Cluster plots...")
all_images <- levels(factor(visium_integrated@meta.data[["orig.ident"]]))
cols <- c("#DDCC77", "#E6AB02", "#7CAE00", "#F8766D", "#00C19A", "#6699CC", "#999999", "#E41A1C", "#377EB8")

png(file.path(plot_dir, "SpatialClusterPlots_Visium_v5_CLR.png"), width = 1200, height = 600, units = 'mm', res = 400)
p_spatial_list <- SpatialDimPlot(visium_integrated, images = all_images, crop = FALSE, pt.size.factor = 1.2, combine = FALSE, cols = cols)
p_spatial_grid <- cowplot::plot_grid(plotlist = p_spatial_list, ncol = 6, nrow = 4) 
print(p_spatial_grid)
dev.off()

# 8. 분석 완료 데이터 저장
message("8. Saving final clustered object...")
saveRDS(visium_integrated, file = file.path(rds_dir, "Final_IntegratedVisium.rds"), compress = TRUE)

cat("\n[성공!] 드디어 모든 파이프라인이 에러 없이 끝났습니다! plots 폴더의 이미지를 확인해주세요.\n")