#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############################### III. INTEGRATION (FINAL FIX) ###################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Fixes included:       Memory optimization, S4/List coercion fix, 
#                       dimnames missing fix (JoinLayers)

library(Seurat)
library(dplyr)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
setwd(base_dir)
rds_dir  <- file.path(base_dir, "rds_OUT")

# 1. 데이터 로드
message("1. Loading preprocessed Visium data...")
stdata_visium <- readRDS(file.path(rds_dir, "SampleListVisium_LogNorm.rds"))

# [데이터 구조 검증]
if (!is.list(stdata_visium)) {
  stop("에러: 데이터가 리스트 형태가 아닙니다. Part 1을 다시 확인해주세요.")
}
if (!all(sapply(stdata_visium, inherits, "Seurat"))) {
  stop("에러: 리스트 내부의 객체가 Seurat 객체가 아닙니다.")
}

# 2. 통합 유전자 선택 (메모리 및 교집합을 위해 2000개로 최적화)
message("2. Selecting integration features...")
features_visium <- SelectIntegrationFeatures(object.list = stdata_visium, nfeatures = 2000)

# 3. 개별 샘플 스케일링 및 PCA 수행
message("3. Running ScaleData & PCA on each sample...")
stdata_visium <- lapply(X = stdata_visium, FUN = function(x) {
  # 경고 방지를 위해 명시적으로 features 전달
  x <- ScaleData(x, features = features_visium, verbose = FALSE)
  x <- RunPCA(x, features = features_visium, verbose = FALSE)
  return(x)
})

# 4. 앵커(Anchors) 찾기
message("4. Finding integration anchors (RPCA)...")
anchors_visium <- FindIntegrationAnchors(
  object.list = stdata_visium,
  normalization.method = "LogNormalize",
  anchor.features = features_visium,
  reduction = "rpca", 
  dims = 1:30
)

# [앵커 생성 검증]
if (is.null(anchors_visium) || !inherits(anchors_visium, "AnchorSet")) {
  stop("에러: 앵커 생성 실패! 샘플 간 유전자 데이터의 교집합이 부족합니다.")
}

# 5. 데이터 통합
message("5. Integrating datasets...")
visium_integrated <- IntegrateData(
  anchorset = anchors_visium,
  normalization.method = "LogNormalize",
  dims = 1:30
)

# [Seurat v5 dimnames 유실 에러 해결]
message("6. Fixing assay structure and layers...")
DefaultAssay(visium_integrated) <- "integrated"

# 데이터가 파편화되어 dimnames가 없어진 경우 JoinLayers로 강제 복구
if (is.null(dimnames(visium_integrated@assays$integrated@data))) {
  message("   -> Joining layers to fix missing dimnames...")
  visium_integrated <- JoinLayers(visium_integrated)
}

# 7. 최종 결과 저장
message("7. Saving integrated data...")
saveRDS(visium_integrated, file = file.path(rds_dir, "IntegratedVisium_LogNorm.rds"), compress = TRUE)

message("\n[축하합니다!] 파트 2 통합이 완벽하게 끝났습니다! 이제 Part 3으로 넘어가세요.")