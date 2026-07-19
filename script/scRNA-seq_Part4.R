#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############################## Cluster Annotation ##############################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted for:          Seurat v5 & Local Windows Environment
# Description:          Load integrated dataset and apply cell type annotations.
#                       Interactive loop removed for seamless pipeline execution.

# Loads required packages
library(Seurat)
library(dplyr)

# [경로 설정]
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"
rds_dir <- file.path(base_dir, "rds_OUT")

# 1. 3단계 최종 결과물인 Seurat 객체 로드
seurat_path <- file.path(rds_dir, "IntgratedRO_LogNorm_1stLvl.rds")
if (!file.exists(seurat_path)) {
  stop("3단계 결과 파일이 존재하지 않습니다. 먼저 Part3을 실행해주세요: ", seurat_path)
}
message("Loading Integrated Seurat Object...")
scdata.integrated <- readRDS(seurat_path)

# 메인 Assay 설정
DefaultAssay(scdata.integrated) <- "integrated"

# [핵심 수정] RStudio 전체 실행 시 충돌을 일으키는 대화형 탐색 루프(while) 제거
message("\nApplying 1st level cell type annotations...")

# Renames cluster identities to their putative cell types.
# 클러스터 ID(0~9)를 분석된 실제 세포 타입 이름으로 일괄 매핑합니다.
scdata.integrated$first.lvl.annot <- recode(scdata.integrated$seurat_clusters,
                                            "0" = "KC", "1" = "Lymphoid",
                                            "2" = "Myeloid", "3" = "FB",
                                            "4" = "MEL", "5" = "VEC",
                                            "6" = "SMC", "7" = "MC",
                                            "8" = "LEC", "9" = "PDC")

# 주클러스터(Ident)를 방금 정의한 세포 타입 라벨로 변경해 줍니다.
Idents(scdata.integrated) <- "first.lvl.annot"

# Writes the current Seurat object to a .rds file.
output_file <- file.path(rds_dir, "IntgratedRO_LogNorm_1stLvl_annot.rds")
saveRDS(scdata.integrated, file = output_file, compress = TRUE)

cat("\n[완료] 4단계 세포 타입 라벨링(Annotation)이 멈춤 없이 완벽하게 완료되었습니다!\n최종 객체 저장 경로: ", output_file, "\n")