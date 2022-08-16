###Mouse hippocampus###
###Data integration & annotation of entire clusters ###

library(Seurat)
library(devtools)
library(Seurat)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(slingshot)
library(SAVER)
library(multtest)
library(metap)
library(tradeSeq)
library(enrichR)

# Before installing SAVER, please install older version 'glmnet' using install_version("glmnet", version = "3.0")

#control_1month#

con_1m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("con_1m", barcode.names[,1], sep = '_')


colnames(con_1m) = barcode.names$V1
rownames(con_1m) = feature.names$V2

con_1m <- CreateSeuratObject(counts = con_1m, project = "con_1m", min.cells = 3, min.features = 200)

con_1m@meta.data$condition <- "control"
con_1m@meta.data$time <- "con_1M"
con_1m@meta.data$time_2 <- "1M"

con_1m[["percent.mt"]] <- PercentageFeatureSet(con_1m, pattern = "^mt-")
head(con_1m@meta.data, 5)
VlnPlot(con_1m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

con_1m <- subset(con_1m, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

#control_3month#

con_3m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("con_3m", barcode.names[,1], sep = '_')


colnames(con_3m) = barcode.names$V1
rownames(con_3m) = feature.names$V2

con_3m <- CreateSeuratObject(counts = con_3m, project = "con_3m", min.cells = 3, min.features = 200)

con_3m@meta.data$condition <- "control"
con_3m@meta.data$time <- "con_3M"
con_3m@meta.data$time_2 <- "3M"


con_3m[["percent.mt"]] <- PercentageFeatureSet(con_3m, pattern = "^mt-")
head(con_3m@meta.data, 5)
VlnPlot(con_3m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

con_3m <- subset(con_3m, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

#control_6month#

con_6m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("con_6m", barcode.names[,1], sep = '_')

colnames(con_6m) = barcode.names$V1
rownames(con_6m) = feature.names$V2

con_6m <- CreateSeuratObject(counts = con_6m, project = "con_6m", min.cells = 3, min.features = 200)

con_6m@meta.data$condition <- "control"
con_6m@meta.data$time <- "con_6M"
con_6m@meta.data$time_2 <- "6M"

con_6m[["percent.mt"]] <- PercentageFeatureSet(con_6m, pattern = "^mt-")
head(con_6m@meta.data, 5)
VlnPlot(con_6m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

con_6m <- subset(con_6m, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)

#app_1month#

app_1m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("app_1m", barcode.names[,1], sep = '_')

colnames(app_1m) = barcode.names$V1
rownames(app_1m) = feature.names$V2

app_1m <- CreateSeuratObject(counts = app_1m, project = "app_1m", min.cells = 3, min.features = 200)

app_1m@meta.data$condition <- "app"
app_1m@meta.data$time <- "app_1M"
app_1m@meta.data$time_2 <- "1M"

app_1m[["percent.mt"]] <- PercentageFeatureSet(app_1m, pattern = "^mt-")
head(app_1m@meta.data, 5)
VlnPlot(app_1m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

app_1m <- subset(app_1m, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 10)


#app_3month#

app_3m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("app_3m", barcode.names[,1], sep = '_')

colnames(app_3m) = barcode.names$V1
rownames(app_3m) = feature.names$V2

app_3m <- CreateSeuratObject(counts = app_3m, project = "app_3m", min.cells = 3, min.features = 200)

app_3m@meta.data$condition <- "app"
app_3m@meta.data$time <- "app_3M"
app_3m@meta.data$time_2 <- "3M"

app_3m[["percent.mt"]] <- PercentageFeatureSet(app_3m, pattern = "^mt-")
head(app_3m@meta.data, 5)
VlnPlot(app_3m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

app_3m <- subset(app_3m, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

#app_6month#

app_6m <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names[,1] <- paste("app_6m", barcode.names[,1], sep = '_')

colnames(app_6m) = barcode.names$V1
rownames(app_6m) = feature.names$V2

app_6m <- CreateSeuratObject(counts = app_6m, project = "app_6m", min.cells = 3, min.features = 200)

app_6m@meta.data$condition <- "app"
app_6m@meta.data$time <- "app_6M"
app_6m@meta.data$time_2 <- "6M"

app_6m[["percent.mt"]] <- PercentageFeatureSet(app_6m, pattern = "^mt-")
head(app_6m@meta.data, 5)
VlnPlot(app_6m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

app_6m <- subset(app_6m, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)

#integration#

app_list <- merge(x = con_1m, y = list(con_3m, con_6m, app_1m, app_3m, app_6m))
app_list <- SplitObject(app_list, split.by = c("time"))
app_list <- lapply(X = app_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = app_list)
app_list <- lapply(X = app_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = app_list, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)
app_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
app_integrated <- ScaleData(app_integrated, verbose = FALSE)
app_integrated <- RunPCA(app_integrated, verbose = FALSE)
app_integrated <- RunUMAP(app_integrated, dims = 1:40)
app_integrated <- FindNeighbors(app_integrated, reduction = "pca", dims = 1:40)
app_integrated <- FindClusters(app_integrated, resolution = 0.3)
DimPlot(app_integrated, reduction = "umap", label = TRUE, label.size = 6)
DimPlot(app_integrated, reduction = "umap", label = TRUE, label.size = 0, split.by = "condition")

#annotation_entire_clusters#

new.cluster.ids <- as.matrix(read.csv("./entire_markers.csv"))
names(new.cluster.ids) <- levels(app_integrated)
app_integrated <- RenameIdents(app_integrated, new.cluster.ids)
levels(app_integrated) <- c("Microglia", "Oligodendrocytes", "Endothelial cells", "ARPs","Astrocytes", "Neural stem cells", "Neuroblasts", "mNEUR", "Pericytes", "Macrophages", "NEUTs", "CPCs","EGs", "Blood cells", "VSMCs", "VLMCs", "Tancytes", "Monocytes", "EC-like cells (Cldn5)", "Unknown glial cells (Ttr, Ecrg4)", "Unknown cells 1 (Cxcl12)", "Unknown cells 2", "Unknown cells 3")
DimPlot(app_integrated, reduction = "umap", pt.size = 0.3, label = F)

