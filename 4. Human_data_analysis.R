###Human data###


c3c9c11a3a12a13

###control_1

C3 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("C3", barcode.names[,1], sep = '_')


colnames(C3) = barcode.names$V1
rownames(C3) = feature.names$V2

C3 <- CreateSeuratObject(counts = C3, project = "C3", min.cells = 3, min.features = 200)

C3@meta.data$condition <- "control"
C3@meta.data$sample <- "C3"


C3[["percent.mt"]] <- PercentageFeatureSet(C3, pattern = "^MT-")
head(C3@meta.data, 5)
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

C3 <- subset(C3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

##control_2

C9 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("C9", barcode.names[,1], sep = '_')


colnames(C9) = barcode.names$V1
rownames(C9) = feature.names$V2

C9 <- CreateSeuratObject(counts = C9, project = "C9", min.cells = 3, min.features = 200)

C9@meta.data$condition <- "control"
C9@meta.data$sample <- "C9"


C9[["percent.mt"]] <- PercentageFeatureSet(C9, pattern = "^MT-")
head(C9@meta.data, 5)
VlnPlot(C9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

C9 <- subset(C9, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

###control_3

C11 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("C11", barcode.names[,1], sep = '_')


colnames(C11) = barcode.names$V1
rownames(C11) = feature.names$V2

C11 <- CreateSeuratObject(counts = C11, project = "C11", min.cells = 3, min.features = 200)

C11@meta.data$condition <- "control"
C11@meta.data$sample <- "C11"


C11[["percent.mt"]] <- PercentageFeatureSet(C11, pattern = "^MT-")
head(C11@meta.data, 5)
VlnPlot(C11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

C11 <- subset(C11, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

###ad_1

AD3 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("AD3", barcode.names[,1], sep = '_')


colnames(AD3) = barcode.names$V1
rownames(AD3) = feature.names$V2

AD3 <- CreateSeuratObject(counts = AD3, project = "AD3", min.cells = 3, min.features = 200)

AD3@meta.data$condition <- "AD"
AD3@meta.data$sample <- "AD3"


AD3[["percent.mt"]] <- PercentageFeatureSet(AD3, pattern = "^MT-")
head(AD3@meta.data, 5)
VlnPlot(AD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AD3 <- subset(AD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

###ad_2

AD12 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("AD12", barcode.names[,1], sep = '_')


colnames(AD12) = barcode.names$V1
rownames(AD12) = feature.names$V2

AD12 <- CreateSeuratObject(counts = AD12, project = "AD12", min.cells = 3, min.features = 200)

AD12@meta.data$condition <- "AD"
AD12@meta.data$sample <- "AD12"


AD12[["percent.mt"]] <- PercentageFeatureSet(AD12, pattern = "^MT-")
head(AD12@meta.data, 5)
VlnPlot(AD12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AD12 <- subset(AD12, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

###ad_3

AD13 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("AD13", barcode.names[,1], sep = '_')


colnames(AD13) = barcode.names$V1
rownames(AD13) = feature.names$V2

AD13 <- CreateSeuratObject(counts = AD13, project = "AD13", min.cells = 3, min.features = 200)

AD13@meta.data$condition <- "AD"
AD13@meta.data$sample <- "AD13"


AD13[["percent.mt"]] <- PercentageFeatureSet(AD13, pattern = "^MT-")
head(AD13@meta.data, 5)
VlnPlot(AD13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AD13 <- subset(AD13, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

##integration##

human_list <- merge(x = C3, y = list(C9,C11,AD3,AD12,AD13))
human_list <- SplitObject(human_list, split.by = "condition")
human_list <- lapply(X = human_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = human_list)
human_list <- lapply(X = human_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = human_list, reference = c(1, 2), reduction = "rpca", dims = 1:50)
AD_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
AD_integrated <- ScaleData(AD_integrated, verbose = FALSE)
AD_integrated <- RunPCA(AD_integrated, verbose = FALSE)
AD_integrated <- RunUMAP(AD_integrated, dims = 1:5)
AD_integrated <- FindNeighbors(AD_integrated, reduction = "pca", dims = 1:5)
AD_integrated <- FindClusters(AD_integrated, resolution = 0.4)

DimPlot(AD_integrated, reduction = "umap", label = TRUE, label.size = 8)

###sort oligodendrocyte###

sort_oligo <- subset(AD_integrated, idents = c("0","1"))
sort_oligo <- ScaleData(sort_oligo, verbose = FALSE)
sort_oligo <- RunPCA(sort_oligo, npcs = 30, verbose = FALSE)
sort_oligo <- RunUMAP(sort_oligo, reduction = "pca", dims = 1:6)
sort_oligo <- FindNeighbors(sort_oligo, reduction = "pca", dims = 1:6)
sort_oligo <- FindClusters(sort_oligo, resolution = 0.7)
DimPlot(sort_oligo, reduction = "umap", split.by = "sample", pt.size = 0.4, label = T, label.size = 8)
DimPlot(sort_oligo, reduction = "umap", split.by = "condition", pt.size = 0.4, label = F, label.size = 8)
DimPlot(sort_oligo, reduction = "umap", label = T, pt.size = 0.9)

VlnPlot(object = sort_oligo, features = c("CD74"), pt.size = 0, combine = T, assay = "RNA", ncol = 1, split.by = "condition")

sort_oligo <- RenameIdents(sort_oligo,'0' = "C",'1' = "D",'2' = "F",'3' = "B",'4' = "C",'5' = "A",'6' = "B",'7' = "G",'8' = "A",'9' = "E")
levels(sort_oligo) <- c("A", "B", "C", "D","E", "F")
DimPlot(sort_oligo, reduction = "umap", label = F, label.size = 7, pt.size = 1.1)

###human_mouse matching###

DAO_gene <- read.csv("./mouse_cell_types_markers.csv")

COP_marker <- DAO_gene$COP
MFOL_marker <- DAO_gene$MFOL
MOL1_marker <- DAO_gene$MOL1
MOL2_marker <- DAO_gene$MOL2
MOL3_marker <- DAO_gene$MOL3
MOL4_marker <- DAO_gene$MOL4
DAO_marker <- DAO_gene$DAO

MFOL_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(MFOL_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
MFOL_marker_predicted$data
write.csv(MFOL_marker_predicted$data, "./MFOL_predicted.csv")

MOL1_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(MOL1_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
MOL1_marker_predicted$data
write.csv(MOL1_marker_predicted$data, "./MOL1_predicted.csv")

MOL2_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(MOL2_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
MOL2_marker_predicted$data
write.csv(MOL2_marker_predicted$data, "./MOL2_predicted.csv")

MOL3_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(MOL3_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
MOL3_marker_predicted$data
write.csv(MOL3_marker_predicted$data, "./MOL3_predicted.csv")

MOL4_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(MOL4_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
MOL4_marker_predicted$data
write.csv(MOL4_marker_predicted$data, "./MOL4_predicted.csv")

DAO_marker_predicted <- DotPlot(sort_oligo, features = as.character(unique(DAO_marker)),cols = c("blue", "red"), dot.scale = 8, assay = "RNA", col.max = 1, col.min = -1) + 
  RotatedAxis()
DAO_marker_predicted$data
write.csv(DAO_marker_predicted$data, "./DAO_predicted.csv")