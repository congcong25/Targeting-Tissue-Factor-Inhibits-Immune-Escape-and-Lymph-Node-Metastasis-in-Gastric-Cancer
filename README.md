# EpiHR_LNM Analysis Scripts

Bulk and single-cell RNA-seq analysis scripts for the EpiHR_LNM project, covering epithelial cells, lymph-node metastasis (LNM), and survival.

## Script Organization

- `main_scripts.R`: The analysis workflow, from meta-cohort screening to the human epithelial and mouse T/NK figures.
- `fun.R`: Helper functions sourced by `main_scripts.R`: `TOP_N()`, `dif_test()`, `surv()`, `ggsur_to_best()`, `ggsur_fun_cag()`, and `gradual_volcano1()`.

## Analysis Workflow

1. Screen the meta-cohort for LNM- and survival-associated genes (`dif_test()`, `surv()`) and save the results as Excel tables.
2. Build the `Gene_all`, `Gene_Epi`, and `Gene_LNM` modules, add the `ALL_HR` and `LNM_HR` scores, and draw the UMAP feature plots.
3. Combine the ssGSEA scores with the clinical data and draw Kaplan-Meier curves.
4. Analyze the human epithelial compartment: Ro/e, `Epi_C1_CCL20` markers, volcano plots, KEGG enrichment, and F3 expression.
5. Analyze the mouse data: cell-type composition, lineage markers, T/NK subsets, correlation with F3, and Treg survival.

## Data and Outputs

The scripts expect the processed input data in the working directory. Figures and tables are written to `2.0_EpiHR_LNM/`, `2.1.3.EpiHR_LNM_scRNAseq/`, and `F3_mouse/`, which are created automatically.

## Running the Scripts

From the repository root (R >= 4.0):

```r
source("main_scripts.R")
```

`main_scripts.R` sources `fun.R` automatically. Required packages include Seurat, SCP, Startrac, GSVA, clusterProfiler, survival, survminer, openxlsx, and related dependencies.

## Notes

For questions about the scripts or required input files, please contact the authors at <nierc@sysucc.org.cn>.
