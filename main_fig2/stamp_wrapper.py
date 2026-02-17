import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import squidpy as sq
import sctm

def stamp_wrapper_fn(counts,coordinates):
    np.Inf = np.inf

    adata = ad.AnnData(counts, obsm={"spatial": coordinates})
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=50)
    sctm.pp.filter_genes(adata, 0.03,  expression_cutoff_99q = 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=6000, flavor="seurat_v3")

    sq.gr.spatial_neighbors(adata)

    n_topics = 20

    # Only hvgs and fit a total of 20 topics
    model = sctm.stamp.STAMP(
        adata[:, adata.var.highly_variable],
        n_topics = n_topics,
    )

    # uses gpu by default to use cpu use device="cpu"
    model.train(device = "cpu")

    topic_prop = model.get_cell_by_topic()
    return(topic_prop, topic_prop.index.tolist())