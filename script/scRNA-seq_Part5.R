#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
########################### SUB-CLUSTER KERATINOCYTES ##########################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Description:          Subset primary clusters, perform high-resolution 
#                       sub-clustering with safe PCA feature intersection and 
#                       robust outlier filtering.

# Loads required packages
library(Seurat)
library(dplyr)
library(cowplot)
library(rcartocolor)
install.packages("tinter")
library(tinter)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
rds_dir <- file.path(base_dir, "rds_OUT")
plot_dir <- file.path(base_dir, "visualization", "plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

input_file <- file.path(rds_dir, "IntgratedRO_LogNorm_1stLvl_annot.rds")
if (!file.exists(input_file)) {
  stop("4단계 결과 파일이 존재하지 않습니다. 먼저 Part4를 실행해주세요: ", input_file)
}

message("Loading annotated 1st level Seurat object...")
scdata.1st.lvl <- readRDS(input_file)

# ----------------- STEP 1: SUBSET KERATINOCYTES (Cluster 0) ----------------- #
message("Subsetting Keratinocytes (Cluster 'KC')...")

# Part 4에서 Annotation한 세포 타입을 기준으로 "KC"만 추출합니다.
# (만약 seurat_clusters == 0 으로 해야 한다면 Ident를 바꿔서 진행합니다)
Idents(scdata.1st.lvl) <- "first.lvl.annot"
if ("KC" %in% levels(scdata.1st.lvl)) {
  cluster0 <- subset(scdata.1st.lvl, idents = "KC")
} else {
  # 만약 라벨링 방식이 달라 0번으로 남았다면 0번 추출
  Idents(scdata.1st.lvl) <- "seurat_clusters"
  cluster0 <- subset(scdata.1st.lvl, idents = "0")
}

# Re-normalize the subset RNA counts
DefaultAssay(cluster0) <- "RNA"
cluster0 <- NormalizeData(cluster0)

# Find variable features in the RNA assay
cluster0 <- FindVariableFeatures(cluster0, selection.method = "vst", nfeatures = 3500)
var_features_rna <- VariableFeatures(cluster0, assay = "RNA")

# [핵심 방어 1] RNA 변동 유전자 중 integrated 어세이에 실제 존재하는 유전자만 교집합으로 추출하여 PCA 에러 방지
DefaultAssay(cluster0) <- "integrated"
integrated_genes <- rownames(cluster0[["integrated"]])
safe_features <- intersect(var_features_rna, integrated_genes)

# Scale and PCA with safe features
message("Scaling & Running PCA for Step 1...")
cluster0 <- ScaleData(cluster0, features = safe_features)
cluster0 <- RunPCA(cluster0, features = safe_features, reduction.name = "pca", verbose = FALSE)

# Sub-clustering
message("Finding sub-neighbors and sub-clusters...")
cluster0 <- FindNeighbors(cluster0, reduction = "pca", dims = 1:20)
cluster0 <- FindClusters(cluster0, resolution = 0.5, graph.name = "integrated_snn")

message("Running sub-level UMAP...")
cluster0 <- RunUMAP(cluster0, reduction = "pca", dims = 1:20, reduction.name = "umap")

saveRDS(cluster0, file = file.path(rds_dir, "KC_v1_LogNorm_2ndLvl_unfiltered.rds"), compress = TRUE)


# -------------------------- STEP 2: FILTER OUTLIERS ------------------------- #
message("Filtering out outlier clusters...")

# [핵심 방어 2] 현재 존재하는 클러스터 번호만 필터링하여 'out of bounds' 에러 원천 차단
current_clusters <- as.numeric(as.character(unique(cluster0$seurat_clusters)))
clst.to.rm <- c(14, 15, 16, 17)
clst.to.keep <- setdiff(current_clusters, clst.to.rm)

Idents(cluster0) <- "seurat_clusters"
cluster0 <- subset(cluster0, idents = as.character(clst.to.keep))
cluster0$seurat_clusters <- droplevels(cluster0$seurat_clusters) # 빈 껍데기 레벨 제거

DefaultAssay(cluster0) <- "RNA"
cluster0 <- NormalizeData(cluster0)
cluster0 <- FindVariableFeatures(cluster0, selection.method = "vst", nfeatures = 2000)
var_features_rna2 <- VariableFeatures(cluster0, assay = "RNA")

# [핵심 방어 3] 필터링 후 다시 한 번 교집합 유전자로 안전하게 PCA 수행
DefaultAssay(cluster0) <- "integrated"
integrated_genes2 <- rownames(cluster0[["integrated"]])
safe_features2 <- intersect(var_features_rna2, integrated_genes2)

message("Re-clustering after outlier removal...")
cluster0 <- ScaleData(cluster0, features = safe_features2)
cluster0 <- RunPCA(cluster0, features = safe_features2, reduction.name = "pca", verbose = FALSE)

cluster0 <- FindNeighbors(cluster0, reduction = "pca", dims = 1:12)
cluster0 <- FindClusters(cluster0, resolution = 0.1, graph.name = "integrated_snn")
cluster0 <- RunUMAP(cluster0, reduction = "pca", dims = 1:12, reduction.name = "umap")

message("Saving filtered sub-cluster UMAP...")
clst.nr <- length(levels(cluster0$seurat_clusters))
hex <- "#88CCEE"
# 색상 부족 에러 방지
clst.cols <- tinter(hex, steps = max(clst.nr, 1), direction = "shades", adjust = 0.1)

png(file.path(plot_dir, "UMAP_KC_Rohjan_LogNorm_2ndLvl_filtered.png"),
    width = 600, height = 200, units = 'mm', res = 400)

p1 <- DimPlot(cluster0, reduction = "umap", group.by = "disease",
              cols = c("#CC6677", "#88CCEE"), shuffle = TRUE)
p2 <- DimPlot(cluster0, reduction = "umap", label = TRUE,
              repel = TRUE, cols = clst.cols, shuffle = TRUE, label.size = 5)

print(plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(1/2, 1/2)))
dev.off()

saveRDS(cluster0, file = file.path(rds_dir, "KC_v1_LogNorm_2ndLvl_filtered.rds"), compress = TRUE)


# ----------------------- STEP 3: CLUSTER ANNOTATION ----------------------- #
message("Assigning 2nd level cell type labels...")

# [핵심 방어 4] 예상치 못한 서브 클러스터 번호가 등장할 경우를 대비한 .default 값 추가
cluster0$second.lvl.annot <- recode(cluster0$seurat_clusters,
                                    "0" = "sbKC", "1" = "bKC",
                                    "2" = "Cycling", "3" = "other KC",
                                    .default = "other KC")

saveRDS(cluster0, file = file.path(rds_dir, "KC_v1_LogNorm_2ndLvl_filtered_annot.rds"), compress = TRUE)


# ------------------------- STEP 4: SUB-CLUSTER PLOTS ------------------------ #
message("Generating final annotated UMAP plot...")

clst.nr_annot <- length(levels(factor(cluster0$second.lvl.annot)))
clst.cols_annot <- tinter(hex, steps = max(clst.nr_annot, 1), direction = "shades", adjust = 0.1)

Idents(cluster0) <- "second.lvl.annot"

png(file.path(plot_dir, "UMAP_KC_Rohjan_LogNorm_2ndLvl_filtered_annot.png"),
    width = 600, height = 200, units = 'mm', res = 400)

p1_annot <- DimPlot(cluster0, reduction = "umap", group.by = "disease",
                    cols = c("#CC6677", "#88CCEE"), shuffle = TRUE)
p2_annot <- DimPlot(cluster0, reduction = "umap", label = TRUE,
                    repel = TRUE, cols = clst.cols_annot, shuffle = TRUE, label.size = 5)

print(plot_grid(p1_annot, p2_annot, align = "h", ncol = 2, rel_widths = c(1/2, 1/2)))
dev.off()

cat("\n[완료] 5단계 Keratinocyte 서브클러스터링 및 라벨링이 에러 없이 끝났습니다!\n")
cat("최종 결과 저장 경로: ", file.path(rds_dir, "KC_v1_LogNorm_2ndLvl_filtered_annot.rds"), "\n")