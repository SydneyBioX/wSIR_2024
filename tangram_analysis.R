# loading packages
library(reticulate)
library(here)
library(anndata)
library(BiocNeighbors)
library(SpatialExperiment)
library(reshape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(philentropy)
library(gghighlight)
library(viridis)

# Tangram wrapper

## set python environment
pythonPath = paste0(here(), "/pythonVENVs/TANGRAM_analysis/.venv/bin/python3.8")
use_python(pythonPath)

source(paste0(here(), "/tangram_evalMetrics_v4.R"))

## Tangram using genes
filePath = '/dski/nobackup/biostat/datasets/SpatialTranscriptomics/'
sp_file = paste0(filePath, 'seqFISH_SpatialMouseAtlas_mouseEmbryos/AnnData/seqFISH_embryo1.h5ad')
sc_file = paste0(filePath, 'seqFISH_SpatialMouseAtlas_mouseEmbryos/AnnData/seqFISH_embryo2.h5ad')
pred_coords_genes <- tangram_run(sp_file, sc_file)
rownames(pred_coords_genes) = colnames(spe_em2)
colnames(pred_coords_genes) = c("X", "Y")
saveRDS(pred_coords_genes, file = paste0(here(), "/data/predCoord_genes.Rds"))

## Tangram using low dimensions
# load all_dfs - a list of low dimension   
# matrices for embryo1 (m1) and embryo2 (m2)
load(paste0(here(), "/exprs_lowdim.Rdata"))
#renaming columns to overcome index error in tangram
colnames(all_dfs$m1_pca) = c(paste0("PCA", 1:ncol(all_dfs$m1_pca))) 
colnames(all_dfs$m2_pca) = c(paste0("PCA", 1:ncol(all_dfs$m2_pca)))

em1_lowDims = all_dfs[c("m1_wsir", "m1_sir", "m1_pca", "m1_lda")]
em2_lowDims = all_dfs[c("m2_wsir", "m2_sir", "m2_pca", "m2_lda")]

filePath = '/dski/nobackup/biostat/datasets/SpatialTranscriptomics/'
# loading seqFISH mouse embryo spe file
spe = readRDS(paste0(filePath, "seqFISH_SpatialMouseAtlas_mouseEmbryos/spe_mouseEmbryo.Rds"))

spe_em1 <- spe[, spe$embryo == "embryo1"]
spCoords_em1 <- matrix(spatialCoords(spe_em1), nrow = ncol(spe_em1), dimnames = list(NULL, c("X", "Y")))
annots_em1 <- data.frame("cellAnnotations" = as.character(spe_em1$celltype_mapped_refined),
                         "X" = spCoords_em1[, 1],
                         "Y" = spCoords_em1[, 2],
                         row.names = colnames(spe_em1))
for (i in names(em1_lowDims)) {
    ad <- AnnData(
        X = em1_lowDims[[i]],
        obs = annots_em1,
        obsm = list(spatial = spCoords_em1))
    write_h5ad(ad, paste0(here(), "/data/seqFISH_mouseEmbryo_", i, ".h5ad"))
}

spe_em2 <- spe[, spe$embryo == "embryo2"]
spCoords_em2 <- matrix(spatialCoords(spe_em2), nrow = ncol(spe_em2), dimnames = list(NULL, c("X", "Y")))
annots_em2 <- data.frame("cellAnnotations" = as.character(spe_em2$celltype_mapped_refined),
                         "X" = spCoords_em2[, 1],
                         "Y" = spCoords_em2[, 2],
                         row.names = colnames(spe_em2))
for (j in names(em2_lowDims)) {
    ad <- AnnData(
        X = em2_lowDims[[j]],
        obs = annots_em2,
        obsm = list(spatial = spCoords_em2))
    write_h5ad(ad, paste0(here(), "/data/seqFISH_mouseEmbryo_", j, ".h5ad"))
}


### wSIR
sp_wsir_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m1_wsir.h5ad')
sc_wsir_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m2_wsir.h5ad')
pred_coords_wsir <- tangram_run(sp_wsir_file, sc_wsir_file)
rownames(pred_coords_wsir) = colnames(spe_em2)
colnames(pred_coords_wsir) = c("X", "Y")
saveRDS(pred_coords_wsir, file = paste0(here(), "/data/predCoord_wsir.Rds"))

### SIR
sp_sir_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m1_sir.h5ad')
sc_sir_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m2_sir.h5ad')
pred_coords_sir <- tangram_run(sp_sir_file, sc_sir_file)
rownames(pred_coords_sir) = colnames(spe_em2)
colnames(pred_coords_sir) = c("X", "Y")
saveRDS(pred_coords_sir, file = paste0(here(), "/data/predCoord_sir.Rds"))

### PCA
sp_pca_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m1_pca.h5ad')
sc_pca_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m2_pca.h5ad')
pred_coords_pca <- tangram_run(sp_pca_file, sc_pca_file)
rownames(pred_coords_pca) = colnames(spe_em2)
colnames(pred_coords_pca) = c("X", "Y")
saveRDS(pred_coords_pca, file = paste0(here(), "/data/predCoord_pca.Rds"))

### LDA
sp_lda_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m1_lda.h5ad')
sc_lda_file = paste0(here(), '/data/seqFISH_mouseEmbryo_m2_lda.h5ad')
pred_coords_lda <- tangram_run(sp_lda_file, sc_lda_file)
rownames(pred_coords_lda) = colnames(spe_em2)
colnames(pred_coords_lda) = c("X", "Y")
saveRDS(pred_coords_lda, file = paste0(here(), "/data/predCoord_lda.Rds"))


## evaluation metrics
true_coords <- spCoords_em2
celltypes = as.character(spe_em2$celltype_mapped_refined)

### Genes
eval_genes = eval_metrics(true_coords, pred_coords_genes, celltypes, neighbor_num = 20,
                          neighHits = TRUE, JSD = TRUE, spearCoef = TRUE)
mean(eval_genes$spearmanCoef$Correlation)
mean(eval_genes$JSdist)

saveRDS(eval_genes, file = paste0(here(), "/data/eval_genes.Rds"))

### wSIR
eval_wsir = eval_metrics(true_coords, pred_coords_wsir, celltypes, neighbor_num = 20,
                         neighHits = TRUE, JSD = TRUE, spearCoef = TRUE)
saveRDS(eval_wsir, file = paste0(here(), "/data/eval_wsir.Rds"))
mean(eval_wsir$spearmanCoef$Correlation)
mean(eval_wsir$JSdist)

### SIR
eval_sir = eval_metrics(true_coords, pred_coords_sir, celltypes, neighbor_num = 20,
                        neighHits = TRUE, JSD = TRUE, spearCoef = TRUE)
saveRDS(eval_sir, file = paste0(here(), "/data/eval_sir.Rds"))
mean(eval_sir$spearmanCoef$Correlation)
mean(eval_sir$JSdist)

### PCA
eval_pca = eval_metrics(true_coords, pred_coords_pca, celltypes, neighbor_num = 20,
                        neighHits = TRUE, JSD = TRUE, spearCoef = TRUE)
saveRDS(eval_pca, file = paste0(here(), "/data/eval_pca.Rds"))
mean(eval_pca$spearmanCoef$Correlation)
mean(eval_pca$JSdist)

### LDA
eval_lda = eval_metrics(true_coords, pred_coords_lda, celltypes, neighbor_num = 20,
                        neighHits = TRUE, JSD = TRUE, spearCoef = TRUE)
saveRDS(eval_lda, file = paste0(here(), "/data/eval_lda.Rds"))
mean(eval_lda$spearmanCoef$Correlation)
mean(eval_lda$JSdist)


## visualisations
jsd_df <- data.frame(Genes = eval_genes$JSdist,
                     wSIR = eval_wsir$JSdist, 
                     SIR = eval_sir$JSdist, 
                     PCA = eval_pca$JSdist, 
                     LDA = eval_lda$JSdist)
jsd_plot = jsd_df %>%
    melt(variable_name = "Method") %>%
    ggplot(aes(x = Method, y = value, fill = Method)) +
    geom_violin() +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    labs(y = "Jensen-Shannon Distance over 20 neighbours") +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          text = element_text(size = 20))

spear_df <- data.frame(Genes = eval_genes$spearmanCoef$Correlation,
                       wSIR = eval_wsir$spearmanCoef$Correlation, 
                       SIR = eval_sir$spearmanCoef$Correlation, 
                       PCA = eval_pca$spearmanCoef$Correlation,
                       LDA = eval_lda$spearmanCoef$Correlation)
spear_plot = spear_df %>%
    melt(variable_name = "Method") %>%
    ggplot(aes(x = Method, y = value, fill = Method)) +
    geom_violin() +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    labs(y = "Spearman's Rank correlation coefficient") +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          text = element_text(size = 20))

jsd_plot + plot_spacer() + spear_plot + plot_layout(widths = c(1, 0.1, 1))

avgHits_df = data.frame(NN_select = eval_wsir$avg_neighHits$NN_select,
                        Genes = eval_genes$avg_neighHits$AvgHits,
                        wSIR = eval_wsir$avg_neighHits$AvgHits, 
                        SIR = eval_sir$avg_neighHits$AvgHits, 
                        PCA = eval_pca$avg_neighHits$AvgHits, 
                        LDA = eval_lda$avg_neighHits$AvgHits)
p1 = avgHits_df %>%
    melt(id.vars = "NN_select", variable_name = "Method") %>%
    ggplot(aes(x = NN_select, y = as.numeric(value), group = Method, color = Method)) +
    geom_line() +
    geom_point() +
    gghighlight() +
    labs(x = "Number of nearest neighbours", 
         y = "Average number of common neighbours") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size = 15))

propHits_df = data.frame(NN_select = eval_wsir$avg_neighHits$NN_select,
                         Genes = eval_genes$avg_neighHits$prop_AvgHits,
                         wSIR = eval_wsir$avg_neighHits$prop_AvgHits, 
                         SIR = eval_sir$avg_neighHits$prop_AvgHits, 
                         PCA = eval_pca$avg_neighHits$prop_AvgHits, 
                         LDA = eval_lda$avg_neighHits$prop_AvgHits)
p2 = propHits_df %>%
    melt(id.vars = "NN_select", variable_name = "Method") %>%
    ggplot(aes(x = NN_select, y = as.numeric(value), group = Method, color = Method)) +
    geom_line() +
    geom_point() +
    gghighlight() +
    labs(x = "Number of nearest neighbours", 
         y = "Average proportion of common neighbours") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size = 15))
p1 + plot_spacer() + p2 + plot_layout(widths = c(1, 0.1, 1))

coord_df = data.frame(celltype = celltypes, 
                      x_orig = as.numeric(true_coords[, 1]),
                      y_orig = as.numeric(true_coords[, 2]),
                      x_pred_genes = as.numeric(pred_coords_genes[, 1]), 
                      y_pred_genes = as.numeric(pred_coords_genes[, 2]),
                      x_pred_wsir = as.numeric(pred_coords_wsir[, 1]), 
                      y_pred_wsir = as.numeric(pred_coords_wsir[, 2]),
                      x_pred_sir = as.numeric(pred_coords_sir[, 1]), 
                      y_pred_sir = as.numeric(pred_coords_sir[, 2]),
                      x_pred_pca = as.numeric(pred_coords_pca[, 1]), 
                      y_pred_pca = as.numeric(pred_coords_pca[, 2]),
                      x_pred_lda = as.numeric(pred_coords_lda[, 1]), 
                      y_pred_lda = as.numeric(pred_coords_lda[, 2]))
colors = c("#635547", "#8EC792", "#9e6762", "#FACB12", "#3F84AA", "#0F4A9C", "#ff891c", "#EF5A9D", "#C594BF", 
           "#DFCDE4", "#139992", "#65A83E", "#8DB5CE", "#005579", "#C9EBFB", "#B51D8D", "#532C8A", "#8870ad",
           "#cc7818", "#FBBE92", "#EF4E22", "#f9decf", "#c9a997", "#C72228", "#f79083", "#F397C0", "#DABE99", 
           "#c19f70", "#354E23", "#C3C388", "#647a4f", "#CDE088", "#f7f79e", "#F6BFCB", "#7F6874", "#989898", 
           "#1A1A1A", "#FFFFFF", "#e6e6e6", "#77441B", "#F90026", "#A10037", "#DA5921", "#E1C239", "#9DD84A")
g1 = ggplot(coord_df, aes(x = x_orig, y = -y_orig)) + 
    geom_point(size = 1, aes(colour = celltype)) +
    scale_color_manual(values = colors) +
    ggtitle("True coordinates - seqFISH Mouse Embryo 2") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))
g2 = ggplot(coord_df, aes(x = x_pred_genes, y = -y_pred_genes)) + 
    geom_point(size = 1, aes(colour = celltype)) +
    scale_color_manual(values = colors) +
    ggtitle("Predicted coordinates - Genes") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))
g3 = ggplot(coord_df, aes(x = x_pred_wsir, y = -y_pred_wsir)) + 
    geom_point(size = 1, aes(colour = celltype)) + 
    scale_color_manual(values = colors) +
    ggtitle("Predicted coordinates - wSIR") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))
g4 = ggplot(coord_df, aes(x = x_pred_sir, y = -y_pred_sir)) + 
    geom_point(size = 1, aes(colour = celltype)) + 
    scale_color_manual(values = colors) +
    ggtitle("Predicted coordinates - SIR") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))
g5 = ggplot(coord_df, aes(x = x_pred_pca, y = -y_pred_pca)) + 
    geom_point(size = 1, aes(colour = celltype)) + 
    scale_color_manual(values = colors) +
    ggtitle("Predicted coordinates - PCA") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))
g6 = ggplot(coord_df, aes(x = x_pred_lda, y = -y_pred_lda)) + 
    geom_point(size = 1, aes(colour = celltype)) + 
    scale_color_manual(values = colors) +
    ggtitle("Predicted coordinates - LDA") + 
    labs(x = "X_coordinate", y = "Y_coordinate") +
    guides(color = guide_legend(title = "Celltypes", override.aes = list(size = 3))) + 
    theme_classic() + 
    theme(text = element_text(size = 12))

g1 + g2 + g3 + g4 + g5 + g6 + plot_layout(guides = 'collect')

