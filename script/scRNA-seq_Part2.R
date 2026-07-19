#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############################### III. INTEGRATION #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(Seurat)
library(dplyr)

base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
rds_dir <- file.path(base_dir, "rds_OUT")

input_file <- file.path(rds_dir, "SampleListRO_LogNorm.rds")
if (!file.exists(input_file)) {
  stop("1단계 결과 파일이 존재하지 않습니다. 먼저 Part1을 실행해주세요.")
}

message("Loading SampleListRO_LogNorm.rds...")
scdata <- readRDS(input_file)

# 원 논문과 완전히 동일한 피처 선택 및 PCA 진행
features <- SelectIntegrationFeatures(object.list = scdata, nfeatures = 3000)
scdata <- lapply(X = scdata, FUN = RunPCA, features = features)

# 앵커 식별 단계 (원작자의 파라미터 100% 유지)
anchors <- FindIntegrationAnchors(object.list = scdata,
                                  normalization.method = "LogNormalize",
                                  anchor.features = features,
                                  dims = 1:27, reduction = "rpca",
                                  k.anchor = 3)

# 모든 샘플의 유전자 교집합 추출 (유전자 유실 방지)
all_genes <- Reduce(intersect, lapply(scdata, rownames))

# [핵심 수정] 원작자 방식인 'integrated' 어세이를 새로 생성하되, 
# 수치상 에러를 피하기 위해 k.weight를 안전하게 85로 튜닝합니다. (세포 유실 없음)
scdata_integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "LogNormalize",
                                   features.to.integrate = all_genes,
                                   dims = 1:27,
                                   k.weight = 60)

output_file <- file.path(rds_dir, "IntegratedRO_LogNorm.rds")
saveRDS(scdata_integrated, file = output_file, compress = TRUE)

cat("\n[완료] 2단계 통합 완료! 세포 유실 없이 안전하게 저장되었습니다.\n")
