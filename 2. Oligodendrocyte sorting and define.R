###OLG sorting & OLGs annotation###

#sorting oligodendrocytes#

app_olig <- subset(app_integrated, idents = c("Oligodendrocytes"))
app_olig <- RunPCA(app_olig, npcs = 20, verbose = FALSE)
app_olig <- RunUMAP(app_olig, reduction = "pca", dims = 1:4)
app_olig <- FindNeighbors(app_olig, reduction = "pca", dims = 1:4)
app_olig <- FindClusters(app_olig, resolution = 0.4)
DimPlot(app_olig, reduction = "umap", split.by = "time", pt.size = 0.6, label = T)
DimPlot(app_olig, reduction = "umap", label = T, label.size = 7)
plot_grid(p1, p2)

#Running SAVER#
con_1m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="con_1M",])))
con_3m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="con_3M",])))
con_6m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="con_6M",])))
app_1m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="app_1M",])))
app_3m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="app_3M",])))
app_6m_for_saver = subset(app_olig, cells = as.vector(rownames(app_olig@meta.data[app_olig@meta.data$time=="app_6M",])))

con_1m_for_saver <- as(as.matrix(con_1m_for_saver@assays$RNA@data), 'sparseMatrix')
con_1m_saver <- saver(con_1m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

con_3m_for_saver <- as(as.matrix(con_3m_for_saver@assays$RNA@data), 'sparseMatrix')
con_3m_saver <- saver(con_3m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

con_6m_for_saver <- as(as.matrix(con_6m_for_saver@assays$RNA@data), 'sparseMatrix')
con_6m_saver <- saver(con_6m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

app_1m_for_saver <- as(as.matrix(app_1m_for_saver@assays$RNA@data), 'sparseMatrix')
app_1m_saver <- saver(app_1m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

app_3m_for_saver <- as(as.matrix(app_3m_for_saver@assays$RNA@data), 'sparseMatrix')
app_3m_saver <- saver(app_3m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

app_6m_for_saver <- as(as.matrix(app_6m_for_saver@assays$RNA@data), 'sparseMatrix')
app_6m_saver <- saver(app_6m_for_saver, ncores = 48, do.fast = T, estimates.only = TRUE)

con_1m_oligo <- CreateSeuratObject(counts = con_1m_saver, project = "con_1m", min.cells = 3, min.features = 200)
con_1m_oligo@meta.data$condition <- "control"
con_1m_oligo@meta.data$time <- "con_1M"
con_1m_oligo@meta.data$time_2 <- "1M"


con_3m_oligo <- CreateSeuratObject(counts = con_3m_saver, project = "con_3m", min.cells = 3, min.features = 200)
con_3m_oligo@meta.data$condition <- "control"
con_3m_oligo@meta.data$time <- "con_3M"
con_3m_oligo@meta.data$time_2 <- "3M"


con_6m_oligo <- CreateSeuratObject(counts = con_6m_saver, project = "con_6m", min.cells = 3, min.features = 200)
con_6m_oligo@meta.data$condition <- "control"
con_6m_oligo@meta.data$time <- "con_6M"
con_6m_oligo@meta.data$time_2 <- "6M"


app_1m_oligo <- CreateSeuratObject(counts = app_1m_saver, project = "app_1m", min.cells = 3, min.features = 200)
app_1m_oligo@meta.data$condition <- "APP"
app_1m_oligo@meta.data$time <- "app_1M"
app_1m_oligo@meta.data$time_2 <- "1M"


app_3m_oligo <- CreateSeuratObject(counts = app_3m_saver, project = "app_3m", min.cells = 3, min.features = 200)
app_3m_oligo@meta.data$condition <- "APP"
app_3m_oligo@meta.data$time <- "app_3M"
app_3m_oligo@meta.data$time_2 <- "3M"


app_6m_oligo <- CreateSeuratObject(counts = app_6m_saver, project = "app_6m", min.cells = 3, min.features = 200)
app_6m_oligo@meta.data$condition <- "APP"
app_6m_oligo@meta.data$time <- "app_6M"
app_6m_oligo@meta.data$time_2 <- "6M"

app_list <- merge(x = con_1m_oligo, y = list(con_3m_oligo, con_6m_oligo, app_1m_oligo, app_3m_oligo, app_6m_oligo))
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
anchors <- FindIntegrationAnchors(object.list = app_list, reference = c(1, 2), reduction = "rpca", dims = 1:50)
app_olig_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
app_olig_integrated <- ScaleData(app_olig_integrated, verbose = FALSE)
app_olig_integrated <- RunPCA(app_olig_integrated, verbose = FALSE)

app_olig_integrated <- RunUMAP(app_olig_integrated, dims = 1:6)
app_olig_integrated <- FindNeighbors(app_olig_integrated, reduction = "pca", dims = 1:6)
app_olig_integrated <- FindClusters(app_olig_integrated, resolution = 0.45)
DimPlot(app_olig_integrated, reduction = "umap", label = TRUE, label.size = 6)
#remove monocytes expressing C1qb;#
app_olig_integrated <- subset(app_olig_integrated, idents = c("9","7"), invert = TRUE)
