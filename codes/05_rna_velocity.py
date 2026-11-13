#!/usr/bin/env python3
"""RNA fractions and ECM-restricted stromal RNA velocity."""
from pathlib import Path
import numpy as np
import pandas as pd
import loompy

ROOT = Path(__file__).resolve().parent.parent
DATA, OUT = ROOT / "data", ROOT / "results"
OUT.mkdir(exist_ok=True)
metadata = pd.read_csv(DATA / "mouse/cell_metadata.csv", index_col="cell")
suffix = {"S13WT_out": "_4", "S19WT_out": "_5", "S31WC_out": "_6",
          "S34KO_out": "_1", "S35KO_out": "_2", "S7KO_out": "_3"}


def cell_id(value):
    sample, barcode = str(value).split(":")
    return barcode.replace("x", "-1") + suffix[sample]


with loompy.connect(str(DATA / "mouse/merged.loom"), mode="r") as loom:
    names = pd.Index([cell_id(s) for s in loom.ca.CellID])
    total = np.empty((len(names), 3), dtype=float)
    for _, indices, view in loom.scan(axis=1, batch_size=256):
        for j, layer in enumerate(("spliced", "unspliced", "ambiguous")):
            total[indices, j] = np.asarray(view.layers[layer][:, :]).sum(axis=0)
    fractions = total / np.maximum(total.sum(axis=1, keepdims=True), 1)
    table = pd.DataFrame(fractions, index=names, columns=["spliced", "unspliced", "ambiguous"])
    table = table.join(metadata[["orig.ident", "Condition", "celltype_fine"]], how="inner")
    table.to_csv(OUT / "Supp_RNA_fractions_cells.csv", index_label="cell")
    table.groupby(["celltype_fine", "Condition"], observed=True)[["spliced", "unspliced", "ambiguous"]].mean().to_csv(
        OUT / "Supp_RNA_fractions_by_celltype.csv")

import anndata as ad
import scanpy as sc
import scvelo as scv

np.random.seed(123)
adata = ad.read_loom(DATA / "mouse/merged.loom", sparse=True)
adata.obs_names = [cell_id(s) for s in adata.obs_names]
cells = metadata.index[metadata.celltype_fine.isin(["Peri", "Myo", "Fibro"])].intersection(adata.obs_names)
adata = adata[cells].copy()
adata.obs = metadata.loc[cells].copy()
adata.var_names_make_unique()
genes = (DATA / "gene_sets/ECM_mouse.txt").read_text().splitlines()
adata = adata[:, adata.var_names.isin(genes)].copy()
adata.raw = adata.copy()
scv.pp.filter_and_normalize(adata)
sc.pp.pca(adata, n_comps=min(30, min(adata.shape) - 1), random_state=123)
sc.pp.neighbors(adata, random_state=123)
sc.tl.umap(adata, random_state=123)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=4)
scv.settings.figdir = str(OUT)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="celltype_fine", density=3,
                               show=False, save="_stromal.pdf")
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)
adata.write_h5ad(OUT / "stromal_velocity.h5ad")
