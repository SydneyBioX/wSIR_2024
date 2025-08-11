# mouse gastrulations
celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Definitive endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     "Gut tube" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Splanchnic mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Lateral plate mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     "Mixed mesenchymal mesoderm" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "Erythroid" = "#f79083",
                     "Blood progenitors" = "#f9decf",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A",
                     
                     "Unknown" = "#FFFFFF",
                     "Low quality" = "#e6e6e6",
                     
                     # somitic and paraxial types
                     # colour from T chimera paper Guibentif et al Developmental Cell 2021
                     "Cranial mesoderm" = "#77441B",
                     "Anterior somitic tissues" = "#F90026",
                     "Sclerotome" = "#A10037",
                     "Dermomyotome" = "#DA5921",
                     "Posterior somitic tissues" = "#E1C239",
                     "Presomitic mesoderm" = "#9DD84A"
)


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
library(MouseGastrulationData)
library(ggpubr)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
# source("plotReducedDimGirafe.R")
# source("plotReducedDimGirafe3.R")

set.seed(2024)

fig_dir = "/Users/sgha9047/Dropbox/wSIR_project/figures/figure_3/raw/"

subset_lower.tri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}

spe = LohoffSeqFISHData(type = "observed", samples = 1:6)

spe

rownames(spe) <- rowData(spe)$SYMBOL

# remove Low quality cells
spe <- spe[,!spe$celltype %in% "Low quality"]

# remove Xist from gene set
spe <- spe[!rownames(spe) %in% "Xist",]

reducedDim(spe, "spatial") <- as.matrix(spatialCoords(spe)[,c("x", "y")]) %*% matrix(c(1,0,0,-1), 2, 2)
spe <- addPerCellQCMetrics(spe, assay.type = "counts")

# keep only cells with 10 counts in 5 genes
spe <- spe[, spe$total >= 10 & spe$detected >= 5]

# keep only genes expressed in at least 1 percent of cells
spe <- addPerFeatureQCMetrics(spe, assay.type = "counts")
spe <- spe[rowData(spe)$detected > 1,]

sort(table(spe$celltype))

spe <- logNormCounts(spe, exprs_values = "counts")

spe

spe <- runPCA(spe)
spe <- runUMAP(spe)

exprs_raw = t(assay(spe, "logcounts"))
coords_raw = as.data.frame(reducedDim(spe, "spatial"))
samples_raw = factor(spe$embryo)
celltypes_raw = factor(spe$celltype)

start = Sys.time()
out = wSIR(X = exprs_raw,
           coords = coords_raw,
           samples = samples_raw,
           # group = celltypes_raw,
           slices = 30,
           alpha = 4,
           maxDirections = 50,
           varThreshold = 0.99,
           optim_params = FALSE)
end = Sys.time()
end - start

getVarImportance = function(a) diag(a %*% solve(crossprod(a)) %*% t(a))

a <- out$directions
varimportance <- getVarImportance(a)
sort(varimportance)

b <- attributes(reducedDim(spe, "PCA"))$rotation
varimportance_pca <- getVarImportance(b)
sort(varimportance_pca)

reducedDim(spe, "wSIR") <- out[[1]]

spe <- runUMAP(spe, dimred = "wSIR", name = "UMAP_wSIR", BPPARAM = MulticoreParam(workers = 4))

pts = 0

np = NNGraphParam(cluster.fun = "louvain", k = 50)
# np = NNGraphParam()

m = clusterRows(reducedDim(spe, "wSIR"), BLUSPARAM = np)
table(m)
spe$m <- m
plotReducedDim(spe, "spatial", colour_by = "m", point_size = pts, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  ggtitle("Cells coloured by wSIR clusters") + 
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo)

mp = clusterRows(reducedDim(spe, "PCA"), BLUSPARAM = np)
table(mp)
spe$mp <- mp
plotReducedDim(spe, "spatial", colour_by = "mp", point_size = pts, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  ggtitle("Cells coloured by PCA clusters") +
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo)


# calculate BCDC per embryo and per cell type for each method

df_split_wSIR = split.data.frame(reducedDim(spe, "wSIR"), list(spe$embryo, spe$celltype))

df_split_PCA = split.data.frame(reducedDim(spe, "PCA"), list(spe$embryo, spe$celltype))

df_split_spatial = split.data.frame(reducedDim(spe, "spatial"), list(spe$embryo, spe$celltype))

bcdc_wSIR = mapply(Rfast::bcdcor, df_split_wSIR, df_split_spatial)
bcdc_PCA = mapply(Rfast::bcdcor, df_split_PCA, df_split_spatial)

identical(names(bcdc_wSIR), names(bcdc_PCA))

bcdc_wSIR_df = na.omit(data.frame(bcdc_wSIR, bcdc_PCA, t(do.call(cbind, strsplit(names(bcdc_wSIR), "\\.")))))
bcdc_wSIR_df$X2 <- factor(bcdc_wSIR_df$X2, levels = names(sort(tapply(bcdc_wSIR_df$bcdc_wSIR, bcdc_wSIR_df$X2, mean), decreasing = TRUE)))



########## downstream analysis

# bcdc plot
g = ggplot(bcdc_wSIR_df, aes(x = X2)) + 
  geom_point(aes(y = bcdc_PCA), colour = "#FFD39B", size = 3, alpha = 0.6) +
  geom_point(aes(y = bcdc_wSIR), colour = "#0000FF", size = 3, alpha = 0.6) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) + 
  ylab("BCDC") + 
  xlab("")

g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_BCDC_celltypes.pdf"), height = 4, width = 10)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_BCDC_celltypes.png"), height = 4, width = 10)


# plot of clusters

g = plotReducedDim(spe, "UMAP_wSIR", colour_by = "m", point_size = pts, scattermore = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = "", ncol = 2)) +
  # ggtitle("Cells coloured by wSIR clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none")
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_clusters.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_clusters.png"), height = 6, width = 8)

g_leg = as.ggplot(get_legend(g + theme(legend.position = "right")))
g_leg

ggsave(g_leg, file = paste0(fig_dir, "MGA_wSIR_UMAP_clusters_leg.pdf"), height = 6, width = 3)
ggsave(g_leg, file = paste0(fig_dir, "MGA_wSIR_UMAP_clusters_leg.png"), height = 6, width = 3)



g = plotReducedDim(spe, "UMAP_wSIR", colour_by = "celltype", point_size = pts, scattermore = TRUE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) + 
  scale_colour_manual(values = celltype_colours) +
  theme(legend.position = "none") + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes.png"), height = 6, width = 8)

g_leg = as.ggplot(get_legend(g + theme(legend.position = "bottom") + guides(colour = guide_legend(title = "", nrow = 3, override.aes = list(size = 10)))))
g_leg

ggsave(g_leg, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes_leg.pdf"), height = 3, width = 16)
ggsave(g_leg, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes_leg.png"), height = 3, width = 16)


g = plotReducedDim(spe, "UMAP_wSIR", colour_by = "celltype", point_size = pts, scattermore = TRUE,
                   other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) + 
  facet_wrap(~embryo) +
  scale_colour_manual(values = celltype_colours) +
  theme(legend.position = "none") + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes_perEmbryo.pdf"), height = 6, width = 16)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_celltypes_perEmbryo.png"), height = 6, width = 16)





g =   plotReducedDim(spe, "UMAP", colour_by = "mp", point_size = pts, scattermore = TRUE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by PCA clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  theme(legend.position = "none") + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "MGA_PCA_UMAP_clusters.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_PCA_UMAP_clusters.png"), height = 6, width = 8)

g_leg = as.ggplot(get_legend(g + theme(legend.position = "right")))
g_leg

ggsave(g_leg, file = paste0(fig_dir, "MGA_PCA_UMAP_clusters_leg.pdf"), height = 6, width = 8)
ggsave(g_leg, file = paste0(fig_dir, "MGA_PCA_UMAP_clusters_leg.png"), height = 6, width = 8)



g =   plotReducedDim(spe, "UMAP", colour_by = "celltype", point_size = pts, scattermore = TRUE) + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by PCA clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  scale_colour_manual(values = celltype_colours) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none") + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "MGA_PCA_UMAP_celltypes.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_PCA_UMAP_celltypes.png"), height = 6, width = 8)


g =   plotReducedDim(spe, "spatial", colour_by = "m", point_size = pts, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = "")) + 
  # ggtitle("Cells coloured by wSIR clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo) + 
  theme(legend.position = "none")
g


ggsave(g, file = paste0(fig_dir, "MGA_wSIR_spatial_clusters.pdf"), height = 8, width = 10)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_spatial_clusters.png"), height = 8, width = 10)


g =   plotReducedDim(spe, "spatial", colour_by = "mp", point_size = pts, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))+ 
  # ggtitle("Cells coloured by PCA clusters") + 
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + 
  ylab("") + 
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo) + 
  theme(legend.position = "none")
g

ggsave(g, file = paste0(fig_dir, "MGA_PCA_spatial_clusters.pdf"), height = 8, width = 10)
ggsave(g, file = paste0(fig_dir, "MGA_PCA_spatial_clusters.png"), height = 8, width = 10)


# heatmap of clusters
p = pheatmap(log10(unclass(table(spe$m, spe$mp))+1),
             # main = "Number of stromal cells (log10)")
             color = colorRampPalette(c("white", "orange", "red", "black"))(100)
)

g = as.ggplot(p)

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_PCA_clusters_confusion.pdf"), height = 4, width = 6)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_PCA_clusters_confusion.png"), height = 4, width = 6)

# calculate ARI
aricode::ARI(spe$m, spe$mp)
aricode::ARI(spe$m, spe$celltype)
aricode::ARI(spe$mp, spe$celltype)

# spatial plots of the most relevant genes
topGenes = names(sort(varimportance, decreasing = TRUE)[1:50])
topGenes_pca = names(sort(varimportance_pca, decreasing = TRUE)[1:50])

# genes = df[df$Top <= 7 & df$summary.logFC > 0,"gene"]
# length(genes)
gList = sapply(c(topGenes[1:5], topGenes_pca[1:5]), function(gene){
  plotReducedDim(spe, "spatial", colour_by = gene, point_size = pts, other_fields = "embryo", scattermore = FALSE) +
    theme(axis.text = element_blank()) +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_blank()) +
    xlab("") + 
    ylab("") + 
    theme(aspect.ratio = 1) +
    theme(legend.position = "none") + 
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
    facet_wrap(~embryo)
}, simplify = FALSE)

g = wrap_plots(gList, ncol = 2, byrow = FALSE)
g

ggsave(g, file = paste0(fig_dir, "MGA_topGenes_wSIR_PCA_gene_spatial.pdf"), height = 16, width = 14)
ggsave(g, file = paste0(fig_dir, "MGA_topGenes_wSIR_PCA_gene_spatial.png"), height = 16, width = 14)


# GO testing on the top genes

length(topGenes)
cat(topGenes)

geneList = rownames(spe) # the gene universe

getEgo = function(topGenes, geneList) {
  
  gene.df = bitr(geneList, fromType = "SYMBOL",
                 toType = c("ENTREZID","SYMBOL"),
                 OrgDb = org.Mm.eg.db)
  head(gene.df)
  
  # obtain the Entrez IDs that correspond to significantly DE genes
  topGenesEntrez = gene.df[gene.df[,1] %in% topGenes,2]
  length(topGenesEntrez)
  # obtain all Entrez IDs that are able to be converted (the gene universe)
  allGenesEntrez = gene.df[gene.df[,1] %in% rownames(spe), 2]
  length(allGenesEntrez)
  
  # overrepresentation test
  egoAll = enrichGO(gene          = topGenesEntrez,
                    universe      = allGenesEntrez,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    minGSSize     = 20,
                    maxGSSize     = 200,
                    readable      = TRUE,
                    pool          = TRUE)
  
  egoAll <- as.data.frame(egoAll)
  
  return(egoAll)
}

egoAll = getEgo(topGenes = topGenes, geneList = geneList)

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

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_top50_GO_barplot.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_top50_GO_barplot.png"), height = 6, width = 8)



### compare egoAll and egoAll_pca
egoAll_pca = getEgo(topGenes = topGenes_pca, geneList = geneList)

goterms = intersect(rownames(egoAll), rownames(egoAll_pca))
df_go = data.frame(wsir = egoAll[goterms, "pvalue"],
                   pca = egoAll_pca[goterms, "pvalue"],
                   label = egoAll[goterms, "Description"],
                   type = egoAll[goterms, "ONTOLOGY"])

g = ggplot(df_go, aes( y = -log(wsir), x = -log(pca))) + 
  geom_point(aes(colour = type)) + 
  geom_text_repel(aes(label = label), data = subset(df_go, rank(wsir) <= 4 | rank(pca) <= 4),
                  size = 5) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_colour_manual(values = ontcols, drop = FALSE) +
  xlab( "-log(P-value) for top PCA genes") + 
  ylab( "-log(P-value) for top wSIR genes") + 
  guides(colour = guide_legend(title = "", override.aes = list(size = 10))) + 
  theme(axis.title = element_text(size = 20))
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_PCA_top50_GO_scatterplot.pdf"), height = 6, width = 12)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_PCA_top50_GO_scatterplot.png"), height = 6, width = 12)


sessionInfo()

########################
########################
# now look at the endothelium subgroup in wSIR

spe$tmp = ifelse(spe$celltype == "Endothelium", "Endothelium", "other")

g = plotReducedDim(spe[, order(spe$tmp)], "UMAP_wSIR", colour_by = "tmp", point_size = pts, scattermore = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = "", ncol = 2)) +
  # ggtitle("Cells coloured by wSIR clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none") + 
  scale_colour_manual(values = setNames(c(celltype_colours["Endothelium"], "grey90"), c("Endothelium", "other")))
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_Endothelium.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_Endothelium.png"), height = 6, width = 8)

# subset the endothelium group and re-run wSIR

endo = spe[,spe$celltype == "Endothelium"]

exprs_raw = t(assay(endo, "logcounts"))
coords_raw = as.data.frame(reducedDim(endo, "spatial"))
samples_raw = factor(endo$embryo)

start = Sys.time()
out = wSIR(X = exprs_raw,
           coords = coords_raw,
           samples = samples_raw,
           # group = celltypes_raw,
           slices = 30,
           alpha = 4,
           maxDirections = 50,
           varThreshold = 0.99,
           optim_params = FALSE)
end = Sys.time()
end - start

reducedDim(endo, "wSIR") <- out$scores

endo <- runUMAP(endo, dimred = "wSIR", name = "UMAP_wSIR", BPPARAM = MulticoreParam(workers = 4))

pts = 0

np = NNGraphParam(cluster.fun = "louvain", k = 50)
# np = NNGraphParam()

m = clusterRows(reducedDim(endo, "wSIR"),
                BLUSPARAM = np)
table(m)
endo$endo_cluster <- m
plotReducedDim(endo, "spatial", colour_by = "endo_cluster", point_size = 5, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  ggtitle("Cells coloured by wSIR clusters") + 
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo)


endo <- runPCA(endo)
endo <- runUMAP(endo, dimred = "PCA")


mp = clusterRows(reducedDim(endo, "PCA"), BLUSPARAM = np)
table(mp)
endo$endo_cluster_pca <- mp

plotReducedDim(endo, "spatial", colour_by = "endo_cluster_pca", point_size = 5, other_fields = "embryo") + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  ggtitle("Cells coloured by wSIR clusters") + 
  theme(aspect.ratio = 1) +
  facet_wrap(~embryo)

table(endo$endo_cluster, endo$endo_cluster_pca)
aricode::ARI(endo$endo_cluster, endo$endo_cluster_pca)

sm = scoreMarkers(endo, groups = endo$endo_cluster)
# pull out the top5 per group - subset to positive LFC and sort by min.auc

numgenes = 15

res = lapply(sm, function(x) {
  xx = subset(x, mean.logFC.detected > 0)
  xx = xx[order(xx$min.AUC, decreasing = TRUE),]
  rownames(xx)[1:numgenes]
})

# now also take out the markers for endothelium
sm_all = scoreMarkers(spe, groups = spe$celltype)[["Endothelium"]]
sm_all = subset(sm_all, mean.logFC.detected > 0)
sm_all = sm_all[order(sm_all$min.AUC, decreasing = TRUE),]
rownames(sm_all)[1:numgenes]

endomarkers = c(rownames(sm_all)[1:numgenes], unlist(res))
duplicated(endomarkers)
endomarkers


endo$endo_cluster <- factor(endo$endo_cluster, levels = rev(levels(endo$endo_cluster)))

g = plotDots(endo, features = unique(endomarkers), group = "endo_cluster", center = TRUE, scale = TRUE) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ylab("") + 
  xlab("Endothelium wSIR subcluster")
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_Endothelium_subclustered_dotplot.pdf"), height = 5, width = 12)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_Endothelium_subclustereddotplot.png"), height = 5, width = 12)




g = plotReducedDim(endo, "UMAP_wSIR", colour_by = "endo_cluster", point_size = 1, scattermore = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = "", ncol = 2)) +
  # ggtitle("Cells coloured by wSIR clusters") +
  theme(axis.text = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  theme(aspect.ratio = 1) +
  theme(legend.position = "right") + 
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 10), title = "")) +
  # scale_colour_manual(values = setNames(c(celltype_colours["Endothelium"], "grey90"), c("Endothelium", "other"))) + 
  NULL
g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_Endothelium_subclustered.pdf"), height = 6, width = 8)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_UMAP_Endothelium_subclustered.png"), height = 6, width = 8)


# spatial plot of the clusters

spe$endo_cluster <- NA
colData(spe)[colnames(endo), "endo_cluster"] <- endo$endo_cluster
spe$endo_cluster <- factor(spe$endo_cluster)


g = plotReducedDim(spe[, order(spe$endo_cluster, na.last = FALSE)], "spatial", colour_by = "endo_cluster",
               point_size = 1, other_fields = "embryo", size_by = "endo_cluster",
               point_shape = 16) + 
  guides(colour = guide_legend(override.aes = list(size = 10))) + 
  facet_wrap(~embryo) + 
  theme(aspect.ratio = 1) +
  scale_size_manual(na.value = 1e-4, values = setNames(rep(1.5, length(unique(spe$endo_cluster))), unique(spe$endo_cluster))) + 
  xlab("") + 
  ylab("") +
  theme(axis.line = element_blank()) + 
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) + 
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(size = 20)) +
  guides(colour = guide_legend(title = "", override.aes = list(size = 10))) + 
  guides(size = guide_none()) + 
  theme(legend.position = "none")

g

ggsave(g, file = paste0(fig_dir, "MGA_wSIR_spatial_Endothelium_subclustered.pdf"), height = 8, width = 12)
ggsave(g, file = paste0(fig_dir, "MGA_wSIR_spatial_Endothelium_subclustered.png"), height = 8, width = 12)


table(endo$embryo, endo$endo_cluster)
