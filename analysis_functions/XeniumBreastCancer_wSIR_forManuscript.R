# xenium breast cancer

library(SingleCellExperiment)
library(scater)
library(StabMap)
library(zellkonverter)
library(wSIR)
library(patchwork)
library(BiocParallel)
library(bluster)
library(pheatmap)
library(scran)
library(dplyr)
library(ggrepel)
library(ggplotify)

set.seed(2024)

getVarImportance = function(a) diag(a %*% solve(crossprod(a)) %*% t(a))

source("/Users/sgha9047/Dropbox/Backup/Projects/wSIR/wSIR_analysis/scripts/plotReducedDimGirafe3.R")

fig_dir = "/Users/sgha9047/Dropbox/wSIR_project/figures/figure_3/raw/"

subset_lower.tri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}

spe_raw <- readH5AD("../../data/xenium_breast_cancer/xenium_rep1.h5ad")
spe_raw <- addPerCellQCMetrics(spe_raw, assay.type = "X")
spe_raw <- spe_raw[, spe_raw$total >= 10 & spe_raw$detected >= 5]

sort(table(spe_raw$Cluster))
# subset to only one type of cells, stromal is the most abundant

spe_raw <- logNormCounts(spe_raw, exprs_values = "X")

spe <- spe_raw[,spe_raw$Cluster == "Stromal"]
# spe <- spe_raw

spe <- addPerFeatureQCMetrics(spe, assay.type = "X")

spe <- runPCA(spe)
spe <- runUMAP(spe)

exprs_raw = t(assay(spe, "logcounts"))
coords_raw = as.data.frame(reducedDim(spe, "spatial"))

start = Sys.time()
out = wSIR(X = exprs_raw,
           coords = coords_raw,
           slices = 30,
           alpha = 4,
           maxDirections = 30,
           varThreshold = 0.99,
           optim_params = FALSE)
end = Sys.time()
end - start

varimportance = getVarImportance(out$directions)
sort(varimportance)

pheatmap(out$W, cluster_rows = FALSE, cluster_cols = FALSE,
         colorRampPalette(c("white", "red"))(100))

reducedDim(spe, "wSIR") <- out[[1]]

spe <- runUMAP(spe, dimred = "wSIR", name = "UMAP_wSIR", BPPARAM = MulticoreParam(workers = 4))

pts = 0

set.seed(2024)
m = clusterRows(reducedDim(spe, "wSIR"),
                BLUSPARAM = TwoStepParam(first = KmeansParam(centers = 1000, iter.max = 10, nstart = 5),
                                         second = NNGraphParam(cluster.fun = "louvain", k = 20)))
table(m)
spe$m <- m
plotReducedDim(spe, "spatial", colour_by = "m", point_size = pts) + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  ggtitle("Cells coloured by wSIR clusters")

set.seed(2024)
mp = clusterRows(reducedDim(spe, "PCA"),
                 BLUSPARAM = TwoStepParam(first = KmeansParam(centers = 1000, iter.max = 10, nstart = 5),
                                          second = NNGraphParam(cluster.fun = "louvain", k = 20)))
table(mp)
spe$mp <- mp
plotReducedDim(spe, "spatial", colour_by = "mp", point_size = pts) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))+ 
  ggtitle("Cells coloured by PCA clusters")


########## downstream analysis 

pcaimportance = getVarImportance(attr(reducedDim(spe, "PCA"), "rotation"))[rownames(spe)]
sort(pcaimportance)
plot(rowData(spe)$detected, pcaimportance, type = "n")
text(rowData(spe)$detected, pcaimportance, label = names(pcaimportance))

plot(rowData(spe)$detected, varimportance, type = "n")
text(rowData(spe)$detected, varimportance, label = names(varimportance))

fit = loess(varimportance ~ log(rowData(spe)$detected), span = 0.7)
plot(log(rowData(spe)$detected), varimportance)
points(log(sort(rowData(spe)$detected)), fit$fitted[order(rowData(spe)$detected)], type = "l")
text(log(sort(rowData(spe)$detected)), varimportance[order(rowData(spe)$detected)], label = rownames(spe)[order(rowData(spe)$detected)])
sort(fit$residuals)


# find what is the identity of the cell type closest to the stromal cells
# split by the wSIR subtype

library(BiocNeighbors)

knn_out = queryKNN(reducedDim(spe_raw[ , spe_raw$Cluster != "Stromal" & spe_raw$Cluster != "Unlabeled" ], "spatial"),
                   reducedDim(spe, "spatial"),
                   k = 1)
knn_tab = unclass(table(droplevels(spe_raw[ , spe_raw$Cluster != "Stromal" & spe_raw$Cluster != "Unlabeled" ]$Cluster[knn_out$index[,1]]), spe$m))
heatmap(knn_tab, scale = "row")

heatmap(chisq.test(knn_tab)$residuals, scale = "none")


p = ComplexHeatmap::Heatmap(chisq.test(knn_tab)$residuals, 
                            column_title = "wSIR Clusters",
                            row_title = "Cell types",
                            col = colorRampPalette(c("blue", "white", "red"))(100),
                            heatmap_legend_param = list(title = "Enrichment"),
                            
)

g = as.ggplot(p)

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_celltypes_proximal.pdf"), height = 4, width = 6)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_celltypes_proximal.png"), height = 4, width = 6)










g = plotReducedDim(spe, "spatial", point_size = pts, scattermore = TRUE) +
  # ggtitle("Stromal") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1)
g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_points.pdf"), height = 4, width = 6)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_points.png"), height = 4, width = 6)


# plot of clusters


g = plotReducedDim(spe, "UMAP_wSIR", colour_by = "m", point_size = pts, scattermore = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by wSIR clusters") +
  ggtitle("wSIR") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) 

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_UMAP_clusters.pdf"), height = 4, width = 4)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_UMAP_clusters.png"), height = 4, width = 4)


g =   plotReducedDim(spe, "UMAP", colour_by = "mp", point_size = pts, scattermore = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by PCA clusters") +
  ggtitle("PCA") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) 
g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_PCA_UMAP_clusters.pdf"), height = 4, width = 4)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_PCA_UMAP_clusters.png"), height = 4, width = 4)


g =   plotReducedDim(spe, "spatial", colour_by = "m", point_size = pts, scattermore = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by wSIR clusters") +
  ggtitle("wSIR") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  # theme(legend.position = "none") + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_spatial_clusters.pdf"), height = 4, width = 5)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_spatial_clusters.png"), height = 4, width = 5)


g =   plotReducedDim(spe, "spatial", colour_by = "mp", point_size = pts, scattermore = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))+ 
  # ggtitle("Cells coloured by PCA clusters") + 
  ggtitle("PCA") + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  # theme(legend.position = "none") + 
  NULL

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_PCA_spatial_clusters.pdf"), height = 4, width = 5)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_PCA_spatial_clusters.png"), height = 4, width = 5)


# heatmap of clusters
# p = pheatmap(log10(unclass(table(spe$m, spe$mp))+1),
#              # main = "Number of stromal cells (log10)")
#              color = colorRampPalette(c("white", "orange", "red", "black"))(100)
# )

p = ComplexHeatmap::Heatmap(log10(unclass(table(spe$mp, spe$m))+1), 
                            row_title = "PCA Clusters",
                            column_title = "wSIR Clusters",
                            col = colorRampPalette(c("white", "orange", "red", "black"))(100),
                            heatmap_legend_param = list(title = "Number of\ncells\n(log10)")
)

g = as.ggplot(p)

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_PCA_clusters_confusion.pdf"), height = 4, width = 5)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_PCA_clusters_confusion.png"), height = 4, width = 5)

# calculate ARI
aricode::ARI(spe$m, spe$mp)

# differential expression between clusters in wSIR

# clusts = c(2,4)
# clusts = c(2,3)
clusts = c(3,5)
mk3 = findMarkers(spe[,spe$m %in% clusts],
                  groups = droplevels(spe[,spe$m %in% clusts]$m),
                  direction = "any")
DE_df = cbind(mk3[[as.character(clusts)[1]]], gene = rownames(mk3[[as.character(clusts)[1]]]))
g = ggplot(DE_df, aes(x = summary.logFC, y = -log10(p.value+1e-100))) + 
  geom_point(colour = "grey") + 
  geom_text_repel(aes(label = gene), data = DE_df[DE_df$Top <= 10,]) + 
  xlab("Log Fold Change") + 
  ylab("-log(P-Value)") + 
  ggtitle(paste0("wSIR Cluster ", clusts[1],
                 " vs Cluster ", clusts[2])) + 
  theme_classic()
g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_clusters_2_4_volcano.pdf"), height = 4, width = 8)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_clusters_2_4_volcano.png"), height = 4, width = 8)

# expression dotplot

sm = scoreMarkers(spe, groups = spe$m)
# pull out the top5 per group - subset to positive LFC and sort by min.auc

numgenes = 15

res = lapply(sm, function(x) {
  xx = subset(x, mean.logFC.detected > 0)
  xx = xx[order(xx$min.AUC, decreasing = TRUE),]
  rownames(xx)[1:numgenes]
})

# now also take out the markers for stromal
sm_all = scoreMarkers(spe_raw, groups = spe_raw$Cluster)[["Stromal"]]
sm_all = subset(sm_all, mean.logFC.detected > 0)
sm_all = sm_all[order(sm_all$min.AUC, decreasing = TRUE),]
rownames(sm_all)[1:numgenes]

stromalmarkers = c(rownames(sm_all)[1:numgenes], unlist(res))
duplicated(stromalmarkers)
stromalmarkers


# spe$stromal_cluster <- factor(spe$endo_cluster, levels = rev(levels(endo$endo_cluster)))

g = plotDots(spe, features = unique(stromalmarkers), group = "m", center = TRUE, scale = TRUE) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ylab("") + 
  xlab("Stromal wSIR subcluster") + 
  theme(legend.position = "bottom")
g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_subclustered_dotplot.pdf"), height = 5, width = 12)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_wSIR_subclustered_dotplot.png"), height = 5, width = 12)





# spatial plots of the top 4 genes for each side

# genes = df[df$Top <= 7 & df$summary.logFC > 0,"gene"]
library(dplyr)
genesUp = DE_df |>
  as.data.frame() |>
  subset(sign(summary.logFC) == 1) |>
  head(4) |>
  rownames()

genesDown = DE_df |>
  as.data.frame() |>
  subset(sign(summary.logFC) == -1) |>
  head(4) |>
  rownames()

genes = c(genesDown, genesUp)
genes

length(genes)
gList = sapply(genes, function(gene){
  plotReducedDim(spe_raw, "spatial", colour_by = gene, point_size = pts, scattermore = FALSE) +
    theme(axis.text = element_blank()) +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("") + 
    ylab("") + 
    theme(aspect.ratio = 1) +
    theme(legend.position = "none") + 
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5))
}, simplify = FALSE)

g = wrap_plots(gList, ncol = 2, byrow = FALSE)

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_cluster_2_gene_spatial.pdf"), height = 14, width = 6)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_clusters_2_gene_spatial.png"), height = 14, width = 6)

g = wrap_plots(gList, nrow = 2, byrow = TRUE)

g

ggsave(g, file = paste0(fig_dir, "xenium_stromal_cluster_2_gene_spatial_wide.pdf"), height = 6, width = 14)
ggsave(g, file = paste0(fig_dir, "xenium_stromal_clusters_2_gene_spatial_wide.png"), height = 6, width = 14)




# GO testing on the top genes

# topGenes = df[df$Top <= 30 & df$summary.logFC > 0,"gene"]
# 
# topGenes = DE_df |>
#   as.data.frame() |>
#   subset(sign(summary.logFC) == 1) |>
#   head(30) |>
#   rownames()
# 
# length(topGenes)
# cat(topGenes)


sm = scoreMarkers(spe, groups = spe$m)
# pull out the top5 per group - subset to positive LFC and sort by min.auc

numgenes = 30

res = lapply(sm, function(x) {
  xx = subset(x, mean.logFC.detected > 0)
  xx = xx[order(xx$min.AUC, decreasing = TRUE),]
  rownames(xx)[1:numgenes]
})

library(clusterProfiler)
library(org.Hs.eg.db)

geneList = rownames(spe) # the gene universe

gene.df = bitr(geneList, fromType = "SYMBOL",
               toType = c("ENTREZID","SYMBOL"),
               OrgDb = org.Hs.eg.db)
head(gene.df)

gList = list()

for (i in seq_along(res)) { 
  
  # obtain the Entrez IDs that correspond to significantly DE genes
  topGenesEntrez = gene.df[gene.df[,1] %in% res[[i]],2]
  length(topGenesEntrez)
  # obtain all Entrez IDs that are able to be converted (the gene universe)
  allGenesEntrez = gene.df[gene.df[,1] %in% geneList, 2]
  length(allGenesEntrez)
  
  # overrepresentation test
  egoAll = enrichGO(gene          = topGenesEntrez,
                    universe      = allGenesEntrez,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    minGSSize     = 20,
                    maxGSSize     = 200,
                    readable      = TRUE,
                    pool          = TRUE)
  
  egoAll <- as.data.frame(egoAll)
  
  df = egoAll |>
    arrange(pvalue) |>
    head(10)
  
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df$ONTOLOGY <- factor(df$ONTOLOGY, levels = c("BP", "MF", "CC"))
  
  ontcols = setNames(c("grey", "black", "brown"), c("BP", "MF", "CC"))
  
  g = ggplot(df, aes(y = -log(pvalue), x = Description)) + 
    geom_col(aes(fill = ONTOLOGY), show.legend = TRUE) + 
    ggtitle(paste0("Cluster", i)) +
    coord_flip() + 
    ylab("-log(P-value)") + 
    xlab("") + 
    theme_classic() + 
    theme(legend.position = "bottom") + 
    scale_fill_manual(values = ontcols, drop = FALSE) +
    guides(fill = guide_legend(title = ""))
  g
  
  gList[[i]] <- g
  
  ggsave(g, file = paste0(fig_dir, "xenium_stromal_cluster_", i, "_GO_barplot.pdf"), height = 6, width = 6)
  ggsave(g, file = paste0(fig_dir, "xenium_stromal_clusters_", i, "_GO_barplot.png"), height = 6, width = 6)
}


### examine Visium and Xenium integration for wSIR

# now read in the xenium data (again)
xen <- readH5AD("../../data/xenium_breast_cancer/xenium_rep1.h5ad")
xen <- addPerCellQCMetrics(xen, assay.type = "X")
xen <- xen[, xen$total >= 10 & xen$detected >= 5]

sort(table(xen$Cluster))
# subset to only one type of cells, stromal is the most abundant
# xen <- xen[,xen$Cluster == "Stromal"]
xen <- logNormCounts(xen, exprs_values = "X")


library(DropletUtils)

spat = read.csv("../../data/DCIS_Visium/spatial/tissue_positions.csv", header = TRUE, row.names = 1)
spe = read10xCounts("../../data/DCIS_Visium/filtered_feature_bc_matrix")
colnames(spe) <- spe$Barcode
rownames(spe) <- make.unique(rowData(spe)$Symbol)
reducedDim(spe, "spatial") <- as.data.frame(as.matrix(spat[spe$Barcode, c("array_row", "array_col")]) %*% matrix(c(1,0, 0, 1),2, 2))

spe <- logNormCounts(spe)
plotReducedDim(spe, "spatial", colour_by = "SFRP4")

spe <- addPerFeatureQCMetrics(spe)
spe <- spe[rowData(spe)$detected > 5,]

genes = intersect(rownames(spe), rownames(xen))

length(genes)

# filter features

# stats <- modelGeneVar(spe)
# genes = getTopHVGs(stats, n = 2000)

spe <- spe[genes,]

start = Sys.time()
out = calculatewSIR(x = spe,
                    dimred = "spatial",
                    slices = 20,
                    alpha = 4,
                    maxDirections = 30,
                    varThreshold = 0.99,
                    optim_params = FALSE)
end = Sys.time()
end - start

reducedDim(spe, "wSIR") <- out$scores

spe <- runUMAP(spe, dimred = "wSIR", name = "wSIR_UMAP")

spe <- runPCA(spe)
spe <- runUMAP(spe, dimred = "PCA", name = "PCA_UMAP")

spe <- addPerCellQCMetrics(spe)


spe$visium_cluster = clusterRows(reducedDim(spe, "wSIR"), BLUSPARAM = NNGraphParam(cluster.fun = "louvain", k = 20))


g = plotReducedDim(spe, "wSIR_UMAP", colour_by = "visium_cluster") + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(title = "", override.aes = list(size = 6))) + 
  ggtitle("Visium DCIS")

g


ggsave(g, file = paste0(fig_dir, "visium_wSIR_UMAP_clusters.pdf"), height = 4, width = 4)
ggsave(g, file = paste0(fig_dir, "visium_wSIR_UMAP_clusters.png"), height = 4, width = 4)



g = plotReducedDim(spe, "spatial", colour_by = "visium_cluster",point_shape = 16, point_alpha = 1, point_size = 2) + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  # guides(colour = guide_legend(title = "", override.aes = list(size = 7)))
  theme(legend.position = "none") + 
  ggtitle("Visium DCIS")

g

ggsave(g, file = paste0(fig_dir, "visium_wSIR_spatial_clusters.pdf"), height = 5, width = 5)
ggsave(g, file = paste0(fig_dir, "visium_wSIR_spatial_clusters.png"), height = 5, width = 5)


plotReducedDim(spe, "PCA_UMAP", colour_by = "SFRP1") + ggtitle("PCA")

attributes(reducedDim(spe, "wSIR"))

library(Rfast)

bcdcor(reducedDim(spe, "wSIR"), reducedDim(spe, "spatial"))
bcdcor(reducedDim(spe, "PCA"), reducedDim(spe, "spatial"))

# now project the xenium cells onto the Visium space
proj = wSIR::projectWSIR(wsir = out, newdata = t(assay(xen, "logcounts")[genes,]))

reducedDim(xen, "wSIR_proj") <- proj

xen <- runPCA(xen)
xen <- runUMAP(xen, dimred = "PCA", name = "PCA_UMAP")

ind = sample(ncol(xen), size = 10000)
bcdcor(reducedDim(xen, "wSIR_proj")[ind,], reducedDim(xen, "spatial")[ind,])
bcdcor(reducedDim(xen, "PCA")[ind,], reducedDim(xen, "spatial")[ind,])

xen <- runUMAP(xen, dimred = "wSIR_proj", name = "wSIR_proj_UMAP")

plotReducedDim(xen, "wSIR_proj_UMAP", colour_by = "MMP12", point_size = 0.2)
plotReducedDim(xen, "PCA_UMAP", colour_by = "MMP12", point_size = 0.2)
plotReducedDim(xen, "spatial", colour_by = "MMP12", point_size = 0.2)



g = plotReducedDim(xen, "wSIR_proj_UMAP", colour_by = "Cluster") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(title = "", override.aes = list(size = 7))) + 
  ggtitle("Xenium DCIS mapped onto Visium embedding") +
  NULL

g

ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_UMAP_celltypes.pdf"), height = 8, width = 8)
ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_UMAP_celltypes.png"), height = 8, width = 8)


xen$tmp <- ifelse(xen$Cluster == "Macrophages_1", "mac", "Not")
plotReducedDim(xen, "wSIR_proj_UMAP", colour_by = "tmp") +
  plotReducedDim(xen, "PCA_UMAP", colour_by = "tmp")

plotReducedDim(xen, "wSIR_proj_UMAP", colour_by = "Cluster") + ggtitle("Xenium cells projected into wSIR Visium space") +
  plotReducedDim(xen, "PCA_UMAP", colour_by = "Cluster") + ggtitle("Xenium cells PCA") 

# plotReducedDimGirafe3(xen, c("wSIR_proj_UMAP", "PCA_UMAP", "spatial"))

xen$selected <- ifelse(reducedDim(xen, "wSIR_proj_UMAP")[,1] > 0 &
                         reducedDim(xen, "wSIR_proj_UMAP")[,1] < 4 &
                         reducedDim(xen, "wSIR_proj_UMAP")[,2] > 3, "selected", "not")

plotReducedDim(xen, "wSIR_proj_UMAP", colour_by = "selected")

xen$mac_sub <- ifelse(xen$Cluster == "Macrophages_1", "Macrophages_1", "Other")
xen$mac_sub[xen$selected == "selected"] <- "Macrophages_1 subset"

g = plotReducedDim(xen[,rev(order(xen$mac_sub))], "spatial", colour_by = "mac_sub",
                   point_size = 0.5) + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  # guides(colour = guide_legend(title = "", override.aes = list(size = 7)))
  theme(legend.position = "bottom") + 
  # scale_colour_manual(values = setNames(c("black", "grey90"), c("selected", "not")))
  scale_colour_manual(values = setNames(c("black", "grey90", "blue"), c("Macrophages_1 subset", "Other", "Macrophages_1"))) + 
  guides(colour = guide_legend(title = "", override.aes = list(size = 15), theme = theme(legend.text = element_text(size = 15))))

g

ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_spatial_selected.pdf"), height = 8, width = 8)
ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_spatial_selected.png"), height = 8, width = 8)


# cluster only the macrophages
# library(bluster)
# xen_sub = xen[, xen$Cluster == "Macrophages_1"]
# xen_sub$cluster_wSIR_proj <- clusterRows(reducedDim(xen_sub, "wSIR_proj"), BLUSPARAM = NNGraphParam())
# xen_sub$cluster_PCA <- clusterRows(reducedDim(xen_sub, "PCA"), BLUSPARAM = NNGraphParam())

# identified as cluster 3 in wSIR
# xen_sub$tmp <- ifelse(xen_sub$cluster_wSIR_proj %in% c("3", "6", "8"), "selected", "Not")
# plotReducedDim(xen_sub, "wSIR_proj_UMAP", colour_by = "tmp") 

mk = findMarkers(xen,
                 groups = droplevels(factor(xen$selected)),
                 direction = "any")
DE_df = cbind(mk[["selected"]], gene = rownames(mk[["selected"]]))
g = ggplot(DE_df, aes(x = summary.logFC, y = -log10(p.value+1e-100))) + 
  geom_point(colour = "grey") + 
  geom_text_repel(aes(label = gene), data = DE_df[DE_df$Top <= 10 & DE_df$logFC.not > 0,]) + 
  xlab("Log Fold Change") + 
  ylab("-log(P-Value)") + 
  # ggtitle(paste0("wSIR Cluster ", clusts[1],
  #                " vs Cluster ", clusts[2])) + 
  theme_classic() + 
  ggtitle("Macrophages_1 subset marker genes")
g


ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_cluster_mac_volcano.pdf"), height = 4, width = 6)
ggsave(g, file = paste0(fig_dir, "visium_projected_wSIR_cluster_mac_volcano.png"), height = 4, width = 6)


varimportance = getVarImportance(out$directions)
sort(varimportance)
pcaimportance = getVarImportance(attr(reducedDim(spe, "PCA"), "rotation"))
sort(pcaimportance)
plot(pcaimportance[genes], varimportance[genes], type = "n", xlab = "Relevance for PCA", ylab = "Relevance for wSIR")
text(pcaimportance[genes], varimportance[genes], label = names(pcaimportance[genes]))
abline(v = sort(pcaimportance, decreasing = TRUE)[10])
abline(h = sort(varimportance, decreasing = TRUE)[10])

importance_df = data.frame(PCA = pcaimportance[genes], wSIR = varimportance[genes], gene = genes)

g = ggplot(importance_df, aes(x = PCA, y = wSIR)) + 
  theme_classic() +
  geom_point(colour = "grey") + 
  geom_text_repel(aes(label = gene), data = subset(importance_df, PCA > 0.7 | wSIR > 0.3 | PCA+wSIR > 0.6),
                  size = 2) +
  xlab("PCA gene relevance") +
  ylab("wSIR gene relevance") +
  theme(axis.title = element_text(size = 10)) +
  # theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  # xlab("") + 
  # ylab("") + 
  theme(aspect.ratio = 1) +
  NULL

g


ggsave(g, file = paste0(fig_dir, "visium_importance_scatterplot.pdf"), height = 4, width = 4)
ggsave(g, file = paste0(fig_dir, "visium_importance_scatterplot.png"), height = 4, width = 4)



g = wrap_plots(sapply(c("MMP12", "OPRPN", "POSTN", "PTGDS", "SFRP4"), function(gene) {
  plotReducedDim(spe, "spatial", colour_by = gene, point_shape = 16, point_alpha = 1, point_size = 1.35) + 
    guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
    ggtitle(gene) +
    theme(axis.text = element_blank()) +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("") + 
    ylab("") + 
    theme(aspect.ratio = 1) +
    theme(legend.position = "none")
}, simplify = FALSE), nrow = 1)

g


ggsave(g, file = paste0(fig_dir, "visium_important_spatial.pdf"), height = 8, width = 20)
ggsave(g, file = paste0(fig_dir, "visium_important_spatial.png"), height = 8, width = 20)



g = plotReducedDim(spe, "spatial", colour_by = "MMP12", point_shape = 16, point_alpha = 1, point_size = 1.35) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") + 
  scale_colour_viridis_c(breaks = c(1,6), labels = c("low", "high")) +
  guides(colour = guide_colourbar(title = "")) + 
  theme(legend.ticks = element_blank()) + 
  # theme(legend.text = element_blank()) + 
  NULL

library(ggpubr)
g_leg = as_ggplot(get_legend(g))

g_leg

ggsave(g_leg, file = paste0(fig_dir, "viridis_leg.pdf"), height = 0.5, width = 1)
ggsave(g_leg, file = paste0(fig_dir, "viridis_leg.png"), height = 0.5, width = 1)



library(clusterProfiler)
library(org.Hs.eg.db)

geneList = genes # the gene universe

gene.df = bitr(geneList, fromType = "SYMBOL",
               toType = c("ENTREZID","SYMBOL"),
               OrgDb = org.Hs.eg.db)
head(gene.df)

topGenes = names(sort(sort(varimportance, decreasing = TRUE)))[1:50]

# obtain the Entrez IDs that correspond to significantly DE genes
topGenesEntrez = gene.df[gene.df[,1] %in% topGenes,2]
length(topGenesEntrez)
# obtain all Entrez IDs that are able to be converted (the gene universe)
allGenesEntrez = gene.df[gene.df[,1] %in% geneList, 2]
length(allGenesEntrez)

# overrepresentation test
egoAll = enrichGO(gene          = topGenesEntrez,
                  universe      = allGenesEntrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  minGSSize     = 20,
                  maxGSSize     = 200,
                  readable      = TRUE,
                  pool          = TRUE)

egoAll <- as.data.frame(egoAll)

df = egoAll |>
  arrange(pvalue) |>
  head(10)

df$Description <- factor(df$Description, levels = rev(df$Description))
df$ONTOLOGY <- factor(df$ONTOLOGY, levels = c("BP", "MF", "CC"))

ontcols = setNames(c("grey", "black", "brown"), c("BP", "MF", "CC"))

g = ggplot(df, aes(y = -log(pvalue), x = Description)) + 
  geom_col(aes(fill = ONTOLOGY), show.legend = TRUE) + 
  coord_flip() + 
  ylab("-log(P-value)") + 
  xlab("") + 
  theme_classic() + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = ontcols, drop = FALSE) +
  guides(fill = guide_legend(title = ""))
g

