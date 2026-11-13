# Mouse kidney VDR analysis

Core analyses for *Vitamin D Receptor Deficiency Accelerates Organ-Specific Premature Aging with the Kidney as the Earliest Target*.

All scripts and shared functions are in this directory. Run the R scripts **from `MouseKidney/codes`**. Inputs are read from `../data`; figures and tables are written to `../results`.

| Script | Analysis |
| --- | --- |
| `01_bulk_aging.R` | Mouse bulk VST normalization, tissue/age expression trends and correlation; Fig. 1A-D. |
| `02_mouse_analysis.R` | Mouse atlas, expression scores, ECM mixture, cell proportions, Wilcoxon DE/GO, and PT/endothelial subtypes; Fig. 3 and related supplementary panels. |
| `03_stromal_trajectory.R` | Mouse and human stromal markers, CytoTRACE comparisons, Monocle2 trajectories and ECM dynamics; Fig. 4. |
| `04_human_bulk_aging.R` | Human Nephroseq age trends and correlation; requires the original spreadsheet. |
| `05_rna_velocity.py` | Spliced/unspliced/ambiguous RNA fractions and ECM-restricted stromal RNA velocity. |
| `common.R` | Shared data access, plotting, statistical and trajectory functions. |

```bash
conda env create -f environment.yml
conda activate mousekidney
Rscript 01_bulk_aging.R
Rscript 02_mouse_analysis.R
Rscript 03_stromal_trajectory.R
python 05_rna_velocity.py
# After supplying ../data/bulk/VDR_CDKN1A_AGE.xlsx:
Rscript 04_human_bulk_aging.R
```

Use Seurat 5.1.0, Monocle 2.30.0 and igraph 1.5.1 as specified in `environment.yml`. The shared functions use four cores and seed 12135; subtype and velocity steps set their own seeds. Allow approximately 16 GB RAM. Input files larger than 10 MiB are excluded from Git. Supply them separately under `../data` using the filenames in `../data/manifest.tsv`; this repository does not use Git LFS.

The mouse atlas contains 42,337 cells and 19,114 genes. It includes the final annotations from `sce_meta_new.csv`. The filtered stromal inputs contain 675 mouse cells and 1,868 human cells. Stored senescence, ECM and CytoTRACE scores are used for figure reproduction. Whole-atlas GO uses Wilcoxon DE with adjusted P < 0.05 and |log2FC| > 1; stromal GO retains raw P < 0.05 and log2FC > 1 from the source analysis.

The archived mouse stromal `counts` layer contains SCT-transformed values. Mouse Monocle fitting therefore uses RNA UMIs from the full atlas for the same cells/genes and selects a root state enriched for pericytes. The original mouse Monocle object is unavailable, so the refitted coordinates can differ. Human plots retain the supplied trajectory and use states 1/7 for the ECM panel.

Human bulk analysis additionally requires the original `../data/bulk/VDR_CDKN1A_AGE.xlsx` spreadsheet. Experimental images and animal measurements require their original source data.
