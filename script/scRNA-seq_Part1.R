#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### INTEGRATION AND CELL-TYPE ANNOTATION OF REFERENCE SC-RNA-SEQ DATASET #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Written by:           Juno Kim
# Project:              Spatial Transcriptomics Combined with Single-cell RNA-
#                       sequencing Unravels the Complex Inflammatory Cell
#                       Network in Atopic Dermatitis
# Institute:            Swiss Institute of Allergy and Asthma Research (SIAF)
# Adapted for:          Seurat v5 & Local Windows Environment

# Loads required packages
library(dplyr)
library(Seurat)

# [경로 설정] 사용자의 로컬 윈도우 디렉토리 경로
base_dir <- "C:/Users/wendy/OneDrive/Desktop/atopic_dermatitis_visium/GSM5907096"

# Lists the absolute paths to all sample sub-directories
sample_paths <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# rds_OUT 등 결과 폴더는 샘플 분석 대상에서 제외합니다.
sample_paths <- sample_paths[!grepl("rds_OUT", sample_paths)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
######## II. PREPROCESSING: QC, FILTERING, AND NORMALIZATION (LOG-NORM) ########
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Creates a list (scdata) of Seurat objects
scdata <- list()
for(i in seq_along(sample_paths)) {
  
  print(paste0("Creating Seurat object for sample ", i, " out of ",
               length(sample_paths), "..."))
  
  # Reads the data one-by-one and creates a Seurat object
  x <- Read10X(data.dir = sample_paths[[i]]) %>%
    CreateSeuratObject(min.cells = 3, min.features = 500,
                       project = "RohjanReference")
  
  # Adds an identifier to each sample (Extract folder name)
  x$orig.ident <- gsub(".*(?<=/)", "", sample_paths[i], perl = TRUE)
  
  print(paste0("Filtering cells from ", i, "..."))
  
  # Filters cells with > 10% mitochondrial genes, > 5000 or < 100 unique feature counts
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  
  x <- subset(x, subset = nFeature_RNA > 100 & nFeature_RNA < 5000)
  x <- subset(x, subset = nCount_RNA > 500 & nCount_RNA < 10000)
  x <- subset(x, subset = percent.mt < 10)
  
  # [Seurat v5 수정] slot 대신 layer 인자를 사용하여 카운트 데이터를 가져옵니다.
  counts <- GetAssayData(x, assay = "RNA", layer = "counts")
  
  # Deletes mitochondrial and ribosomal genes
  counts <- counts[-(grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                          rownames(counts))), ]
  counts <- counts[-(grep("^MT-", rownames(counts))), ]
  x <- subset(x, features = rownames(counts))
  
  scdata[[length(scdata)+1]] <- x
  rm(x)
}

# [출력 경로 설정] 출력 폴더가 없을 경우 자동 생성
out_dir <- file.path(base_dir, "rds_OUT")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#----------------------SAVE/LOAD filtered Seurat objects-----------------------#
# Saves list of raw Seurat objects 
saveRDS(scdata, file = file.path(out_dir, "SeuratObjectList_preNorm.rds"), compress = TRUE)
#------------------------------------------------------------------------------#

# Remove samples with low cell numbers (< 250) after filtering
i <- 1; j <- 1
for(sample in scdata){
  count <- length(colnames(sample))
  if(count < 250) {
    scdata <- scdata[-i]
    cat("Sample ", j, " has ", count, " cells.\nSample ", j, " removed.\n")
  } else {
    i = i + 1
    cat("Sample ", j, " has ", count, " cells.\n")
  }
  j = j + 1
  rm(sample)
}

# [Seurat v5 수정] Normalizing and scaling each sample by 10'000
scdata <- lapply(X = scdata, FUN = function(x){
  x <- NormalizeData(x, assay = "RNA",
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x))
})

# Writes the list of normalized Seurat objects to a .rds file
output_file <- file.path(out_dir, "SampleListRO_LogNorm.rds")
saveRDS(scdata, file = output_file, compress = TRUE)

cat("\n[완료] 1단계 전처리가 성공적으로 완료되었습니다!\n결과 파일 저장 경로: ", output_file, "\n")