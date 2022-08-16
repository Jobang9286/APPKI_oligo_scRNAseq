###Pseudotime analysis###

app_olig_integrated_1 <- app_olig_integrated

sce_1 <- as.SingleCellExperiment(app_olig_integrated_1, assay = "RNA")

shuffle <- sample(ncol(sce_1))

sce_1 <- slingshot(sce_1, reducedDim = 'UMAP', clusterLabels = sce_1$ident, start.clus = 'COP+NFOL', approx_points = 150)
layout(matrix(c(1, 1, 2, 3), 2))
par(mar = c(4.5, 4, 1, 1))

plot(reducedDims(sce_1)$UMAP, col = cell_colors, pch=16, asp = 1)
lines(SlingshotDataSet(sce_1), lwd=2, col='black')

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(sce_1$condition, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(sce_1$seurat_clusters, hue_pal())

plot(reducedDims(sce_1)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = .5)[cut(sce_1$slingPseudotime_3, breaks = 100)][shuffle])

lines(SlingshotDataSet(sce_1))

########################################################
#Density plot of DAO lineage; time point
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "1M", 3], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "3M", 3], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "6M", 3], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(4, "Set1")[1], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(4, "Set1")[2], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(4, "Set1")[3], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(4, "Set1")[4], alpha = .5))

ks.test(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "1M", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$time_2 == "3M", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$time_2 == "6M", 3], na.rm = T)

########################################################
#Density plot of DAO lineage; cell type 
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "COP+NFOL", 3], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "MFOL", 3], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "DAO", 3], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(8, "Set1")[4], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(8, "Set1")[5], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(8, "Set1")[6], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(4, "Set1")[4], alpha = .5))

ks.test(slingPseudotime(sce_1)[colData(sce_1)$ident == "COP+NFOL", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$ident == "MFOL", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$ident == "DAO", 3], na.rm = T)

########################################################
#Density plot of DAO lineage; condition
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "control", 3], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "APP", 3], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(8, "Set1")[7], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(8, "Set1")[8], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(8, "Set1")[6], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(4, "Set1")[4], alpha = .5))

ks.test(slingPseudotime(sce_1)[colData(sce_1)$condition == "control", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$condition == "APP", 3], na.rm = T)

########################################################
#Density plot of MOL lineage; time point
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "1M", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "3M", 1], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "6M", 1], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(4, "Set1")[1], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(4, "Set1")[2], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(4, "Set1")[3], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(4, "Set1")[4], alpha = .5))

ks.test(slingPseudotime(sce_1)[colData(sce_1)$time_2 == "1M", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$time_2 == "3M", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$time_2 == "6M", 3], na.rm = T)

########################################################
#Density plot of MOL lineage; cell type
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "COP+NFOL", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "MFOL", 1], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "MOL3", 1], na.rm = T),
           a_4 = density(slingPseudotime(sce_1)[colData(sce_1)$ident == "MOL2", 1], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x, ds$a_4$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y, ds$a_4$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(8, "Set1")[4], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(8, "Set1")[5], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(8, "Set1")[6], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(8, "Set1")[7], alpha = .5))


ks.test(slingPseudotime(sce_1)[colData(sce_1)$ident == "COP+NFOL", 1], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$ident == "MFOL", 1], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$ident == "MOL3", 1], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$ident == "MOL2", 1], na.rm = T)

########################################################
#Density plot of MOL lineage; condition
########################################################

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "control", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "APP", 1], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(10, "Set1")[8], alpha = .5))

polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(10, "Set1")[9], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(8, "Set1")[6], alpha = .5))
polygon(c(min(ds$a_4$x), ds$a_4$x, max(ds$a_4$x)), c(0, ds$a_4$y, 0),
        col = alpha(brewer.pal(4, "Set1")[4], alpha = .5))

ks.test(slingPseudotime(sce_1)[colData(sce_1)$condition == "control", 3], na.rm = T,
        slingPseudotime(sce_1)[colData(sce_1)$condition == "APP", 3], na.rm = T)

########################################################
#Gene expression pattern on pseudotime
########################################################

sce <- sce_1

slingsce<-SlingshotDataSet(sce)
pseudotimeED <- slingPseudotime(slingsce, na = FALSE)
cellWeightsED <- slingCurveWeights(slingsce)
conditions <- colData(sce)$condition
counts<-sce@assays@data@listData$counts
set.seed(3)
sce <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T, genes = 1:2500, conditions = as.factor(sce$condition))
mean(rowData(sce)$tradeSeq$converged)

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))

assocRes <- rowData(sce)$assocRes
assocRes

p1 <- plotSmoothers(sce, assays(sce)$counts, gene = c("App"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))
p2 <- plotSmoothers(sce, assays(sce)$counts, gene = c("Cd74"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))
p3 <- plotSmoothers(sce, assays(sce)$counts, gene = c("Apoe"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))
p4 <- plotSmoothers(sce, assays(sce)$counts, gene = c("Mgst3"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))
p5 <- plotSmoothers(sce, assays(sce)$counts, gene = c("Anxa5"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))
p6 <- plotSmoothers(sce, assays(sce)$counts, gene = c("Klk6"), alpha = 0.6, border = T, lwd = 2, curvesCols = c("#F98D86","#02C0C4"), pointCol = c("#000000","purple"))

plot_grid(p1,p2,p3,p4,p5,p6)