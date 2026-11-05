# =============================================================================
# Project: EpiHR_LNM analysis
# Description: Bulk/single-cell RNA-seq analysis of epithelial cells,
#              lymph-node metastasis and survival.
# Requirements: R >= 4.0; see package list below.
# Run: source("CR_github.R") after placing data in working directory.
# =============================================================================

# ---- 0. Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(openxlsx)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(survival)
  library(survminer)
  library(paletteer)
  library(pheatmap)
  library(qs)
  library(Startrac)
  library(SCP)
  library(ggVolcano)
  library(ggpubr)
  library(ggcorrplot)
  library(ggsci)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSVA)
  library(future)
  library(latex2exp)
  library(ggrepel)
})

# ---- 1. Global settings ----------------------------------------------------
set.seed(123)
options(stringsAsFactors = FALSE)

dirs <- c("2.0_EpiHR_LNM",
          "2.1.3.EpiHR_LNM_scRNAseq",
          "F3_mouse",
          "F3_mouse/3.0_T_cell_main")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

source("fun.R")

# ---- 2. Load data ----------------------------------------------------------
load("meta_data.rdata")          # provides rt_total, clin_total
gene <- rownames(rt_total)

# ---- 3. Meta-cohort: LNM-associated genes ----------------------------------
meta_data_dif <- do.call(rbind, lapply(
  gene, FUN = dif_test, clin = clin_total, rt = rt_total,
  group = "lymph_node_meta"
))
write.xlsx(meta_data_dif,
           file = "2.0_EpiHR_LNM/meta_data_dif_LNM.xlsx",
           sheetName = "dif_LNM", row.names = FALSE)


# ---- 4. Meta-cohort: survival-associated genes -----------------------------
meta_data_sur_conti <- do.call(rbind, lapply(
  gene, FUN = surv, rt = rt_total, clin = clin_total,
  method = "continuous", time = "OS_m", event = "OS_status", minprop = 0.1
))
write.xlsx(meta_data_sur_conti,
           file = "2.0_EpiHR_LNM/meta_data_survival.xlsx",
           sheetName = "continuous", row.names = FALSE)

# ---- 5. Load single-cell objects -------------------------------------------
scobj        <- qs::qread("2.1.3.EpiHR_LNM_scRNAseq/scobj.anno_mal_nor.qs")
scobj_epi <- qs::qread("2.1.3.EpiHR_LNM_scRNAseq/scobj.final.qs")  # final epithelial object

scobj$anno_epi <- ifelse(scobj$anno %in% "Epithelial", "Epi", "TME")

# ---- 6. Epi vs TME markers -------------------------------------------------
Idents(scobj) <- "anno_epi"
anno_epi_markers_ours <- FindAllMarkers(scobj, assay = "RNA",
                                        logfc.threshold = 0, verbose = FALSE)
fwrite(anno_epi_markers_ours,
       "2.1.3.EpiHR_LNM_scRNAseq/epi.markers.csv", row.names = FALSE)

# ---- 7. Malignant epithelial cells: primary vs LNM -------------------------
scobj_Mal <- subset(scobj, anno1 %in% "Mal_epi")
Idents(scobj_Mal) <- "type"
anno_epi_markers_ours_lnm <- FindAllMarkers(scobj_epi, assay = "RNA",
                                            logfc.threshold = 0, verbose = FALSE)
fwrite(anno_epi_markers_ours_lnm,
       "2.1.3.EpiHR_LNM_scRNAseq/anno_epi_markers_lnm.csv", row.names = FALSE)

# ---- 8. Define gene modules ------------------------------------------------
Gene_all <- intersect(
  meta_data_dif$gene_name[meta_data_dif$mean_dif > 0 &
                            meta_data_dif$pvalue_wilcox < 0.05],
  meta_data_sur_conti$gene_name[meta_data_sur_conti$HR > 1 &
                                  meta_data_sur_conti$pvalue_HR < 0.05]
) %>% intersect(rownames(scobj))

anno_epi_markers_ours_lnm <- fread(
  "2.1.3.EpiHR_LNM_scRNAseq/anno_epi_markers_lnm.csv", data.table = FALSE)
anno_epi_markers_ours <- fread(
  "2.1.3.EpiHR_LNM_scRNAseq/epi.markers.csv", data.table = FALSE)

Gene_Epi <- intersect(
  Gene_all,
  anno_epi_markers_ours$gene[
    anno_epi_markers_ours$cluster %in% "Epi" &
      anno_epi_markers_ours$avg_log2FC > 0 &
      anno_epi_markers_ours$p_val_adj < 0.05 &
      anno_epi_markers_ours$pct.1 > 0.2 &
      anno_epi_markers_ours$pct.2 < 0.5 &
      anno_epi_markers_ours$pct.1 > anno_epi_markers_ours$pct.2
  ]
)

Gene_LNM <- intersect(
  Gene_Epi,
  anno_epi_markers_ours_lnm$gene[
    anno_epi_markers_ours_lnm$cluster %in% "MT" &
      anno_epi_markers_ours_lnm$avg_log2FC > 0 &
      anno_epi_markers_ours_lnm$p_val_adj < 0.05 &
      anno_epi_markers_ours_lnm$pct.1 > 0.2 &
      anno_epi_markers_ours_lnm$pct.2 < 0.5 &
      anno_epi_markers_ours_lnm$pct.1 > anno_epi_markers_ours_lnm$pct.2
  ]
)

scobj <- AddModuleScore(scobj,
                        features = list(Gene_all, Gene_LNM),
                        name = c("ALL_HR", "LNM_HR"))

# ---- 9. Figure 1B: ALL_HR FeaturePlot --------------------------------------
p_all_hr <- FeatureDimPlot(
  srt = scobj, features = "ALL_HR1",
  reduction = "umap", theme_use = "theme_blank", assay = "RNA",
  show_stat = FALSE, xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 2, label_insitu = TRUE
)
ggsave("2.0_EpiHR_LNM/ALL_HR_FeaturePlot.pdf", p_all_hr,
       width = 20, height = 20)

# ---- 10. Figure 1C: LNM_HR FeaturePlot -------------------------------------
p_lnm_hr <- FeatureDimPlot(
  srt = scobj, features = "LNM_HR2",
  reduction = "umap", theme_use = "theme_blank", assay = "RNA",
  show_stat = FALSE, xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 2, label_insitu = TRUE
)
ggsave("2.0_EpiHR_LNM/LNM_HR_FeaturePlot.pdf", p_lnm_hr,
       width = 20, height = 20)

# ---- 11. ssGSEA scores and survival ----------------------------------------
Gene_list <- list("Gene_all" = Gene_all, "Gene_LNM" = Gene_LNM)
exp1 <- GSVA::gsva(rt_total, Gene_list, method = "ssgsea", parallel.sz = 8)
gene <- rownames(exp1)

clin_total <- cbind(clin_total[, 1:19], t(exp1))
clin_total$F3 <- as.numeric(rt_total["F3", ])

clin_total$Gene_all_median <- ifelse(
  clin_total$Gene_all > median(clin_total$Gene_all), 1, 0)
clin_total$Gene_LNM_median <- ifelse(
  clin_total$Gene_LNM > median(clin_total$Gene_LNM), 1, 0)
clin_total$Gene_F3_median <- ifelse(
  clin_total$F3 > median(clin_total$F3), 1, 0)

# Figure 1D

pdf("2.0_EpiHR_LNM/Gene_ALL_HR.pdf",width = 8, height = 6)
p1 <- ggsur_fun_cag(total = clin_total, gene = "Gene_all_median",
                    color_palette1 = "d3", x_break = 36, conf.int = FALSE,
                    OS_m = "OS_m", OS_status = "OS_status")
print(p1, newpage = FALSE)
dev.off()


# Figure 1E
pdf("2.0_EpiHR_LNM/Gene_Epi_LNM_HR.pdf",width = 8, height = 6)

p11 <- ggsur_fun_cag(total = clin_total, gene = "Gene_LNM_median",
                     color_palette1 = "d3", x_break = 36, conf.int = FALSE,
                     OS_m = "OS_m", OS_status = "OS_status")
print(p11, newpage = FALSE)
dev.off()

# ---- 12. Figure 1F: epithelial cell UMAP -----------------------------------
p_epi_umap <- CellDimPlot(
  srt = scobj_epi, group.by = "epi_cell", show_stat = FALSE, pt.size = 1.5,
  reduction = "umap", theme_use = "theme_blank",
  xlab = "UMAP1", ylab = "UMAP2"
)
ggsave("2.0_EpiHR_LNM/Dim_plot_epi.pdf", p_epi_umap, width = 10, height = 10)

# ---- 13. Ro/e analysis for epithelial cells --------------------------------
in.dat <- scobj_epi@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = FALSE,
                      colname.cluster = "epi_cell",
                      colname.tissue  = "type",
                      method = "chisq",
                      min.rowSum = 0)

plot_ro_e <- function(R_oe, filename, main = "Ro/e", transpose = FALSE,
                      width = 6, height = 8) {
  mat <- if (transpose) t(R_oe) else R_oe
  pheatmap(mat,
           border = TRUE,
           color = colorRampPalette(c("#0044BB", "white", "#AAAAAA", "red"))(100),
           annotation_names_row = FALSE,
           annotation_names_col = TRUE,
           legend = FALSE,
           annotation_legend = TRUE,
           border_color = NA,
           cellwidth = 30, cellheight = 30,
           display_numbers = TRUE,
           number_format = "%.2f",
           number_color = "black",
           fontsize_number = 8,
           main = main,
           fontsize = 10, fontsize_row = 10, fontsize_col = 10,
           angle_col = "45",
           show_rownames = TRUE, show_colnames = TRUE,
           cluster_cols = FALSE, cluster_rows = FALSE,
           filename = filename, width = width, height = height)
}
plot_ro_e(R_oe, "2.0_EpiHR_LNM/Roe.epi1.pdf", transpose = TRUE)

# ---- 14. Focus on Epi_C1_CCL20 ---------------------------------------------
scobj_epi$C1 <- ifelse(scobj_epi$epi_cell %in% "Epi_C1_CCL20", "Yes", "No")
Idents(scobj_epi) <- "C1"
cluster.epi.markers_C1 <- FindAllMarkers(scobj_epi, assay = "RNA",
                                         min.pct = 0.25,
                                         logfc.threshold = 0,
                                         verbose = FALSE)

DEG2 <- cluster.epi.markers_C1[cluster.epi.markers_C1$cluster %in% "Yes", ]
DEG2$dif_pct <- DEG2$pct.1 - DEG2$pct.2
names(DEG2) <- c("pValue", "Log2FC", "pct.1", "pct.2",
                 "FDR", "cluster", "GeneName", "dif_pct")
rownames(DEG2) <- DEG2$GeneName
DEG2$Log2FC[DEG2$Log2FC < -2.5] <- -2.5

data <- add_regulate(DEG2, log2FC_name = "Log2FC",
                     fdr_name = "FDR", log2FC = 1, fdr = 0.05)
data$row <- rownames(data)

top10 <- cluster.epi.markers_C1 %>% TOP_N(10, pct.1 = 0.2)
top10 <- top10[top10$avg_log2FC > 0 & top10$cluster %in% "Yes", ]

# Figure 1F
p_volcano1 <- gradual_volcano1(
  data, x = "log2FoldChange", y = "padj",
  label = "row", label_number = 10,
  output = FALSE, custom_label = Gene_LNM
) + ggsci::scale_color_gsea() + ggsci::scale_fill_gsea()
ggsave("2.0_EpiHR_LNM/Volcano_C1_CCL20_Gene_LNM.pdf", p_volcano1,
       width = 6, height = 6)

# ---- 15. Figure S1A: cell-type UMAP ----------------------------------------
pal <- paletteer_d("ggsci::nrc_npg")[c(1, 3, 2, 5, 4, 6, 7, 9, 10)]
p_anno <- CellDimPlot(
  srt = scobj, group.by = "anno", show_stat = FALSE,
  reduction = "umap", theme_use = "theme_blank",
  xlab = "UMAP1", ylab = "UMAP2", palcolor = pal
)
ggsave("2.0_EpiHR_LNM/anno.umap.pdf", p_anno,
       width = 7.09, height = 5.6, dpi = 600)

# ---- 16. Figure S1B: split by type -----------------------------------------
p_anno_split <- CellDimPlot(
  srt = scobj, group.by = "anno", show_stat = FALSE, split.by = "type",
  reduction = "umap", theme_use = "theme_blank", palcolor = pal
)
ggsave("2.0_EpiHR_LNM/DimPlot.type.anno.pdf", p_anno_split,
       width = 20, height = 10)

# ---- 17. Figure S1E: volcano for C1 markers --------------------------------
p_volcano2 <- gradual_volcano(
  data, x = "log2FoldChange", y = "padj",
  label = "row", label_number = 10,
  output = FALSE, custom_label = top10$gene
) + ggsci::scale_color_gsea() + ggsci::scale_fill_gsea()
ggsave("2.0_EpiHR_LNM/Volcano_C1_CCL20.pdf", p_volcano2,
       width = 6, height = 6)

# ---- 18. Figure S1F: KEGG enrichment for C1 markers ------------------------
DEG2 <- cluster.epi.markers_C1
ids  <- bitr(DEG2$gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
DEG2 <- merge(DEG2, ids, by.x = "gene", by.y = "SYMBOL")
DEG2 <- DEG2[DEG2$avg_log2FC > 0.25, ]
gcSample <- split(DEG2$ENTREZID, DEG2$cluster)
xx <- compareCluster(gcSample, fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05)
p_kegg <- dotplot(xx, showCategory = 10, includeAll = TRUE)
ggsave("2.0_EpiHR_LNM/kegg_c1_ccl20.pdf", p_kegg, width = 6, height = 8)

# ---- 19. Figure S1G: F3 FeaturePlot on epithelial cells --------------------
scobj_epi$F3_Feature <- as.numeric(scobj_epi@assays$RNA["F3", ])
p_f3_epi <- FeatureDimPlot(
  srt = scobj_epi, features = "F3_Feature",
  reduction = "umap", theme_use = "theme_blank", assay = "RNA",
  show_stat = FALSE, xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 2, label_insitu = TRUE
)
ggsave("2.0_EpiHR_LNM/F3_FeaturePlot_epi.pdf", p_f3_epi,
       width = 20, height = 20)

# ---- 20. Figure 4F: CCL20 by F3 in malignant epithelial cells --------------
scobj_epi_mal <- subset(scobj_epi, epi_cell1 %in% "Mal")
my_comparisons <- list(c("pos", "neg"))

p_vln <- VlnPlot(scobj_epi_mal, features = "CCL20",
                 pt.size = 0, group.by = "F3") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", face = "bold", size = 12),
        axis.text.y = element_text(color = "black", face = "bold"),
        axis.title.y = element_text(color = "black", face = "bold", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1.2,
                                    linetype = "solid"),
        panel.spacing = unit(0.12, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),
        legend.position = "none") +
  stat_compare_means(method = "t.test", hide.ns = FALSE,
                     comparisons = my_comparisons,
                     label = "p.signif", bracket.size = 0.8,
                     tip.length = 0, size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_fill_manual(values = c("#FF5744", "#208A42"))
ggsave("2.0_EpiHR_LNM/scobj_epi_mal_CCL20_by_F3.pdf", p_vln,
       width = 6, height = 8)

# ---- 21. Figure S5F: markers of F3 pos vs neg ------------------------------
Idents(scobj_epi_mal) <- "F3"
cluster.markers_F3 <- FindAllMarkers(scobj_epi_mal, assay = "RNA",
                                     logfc.threshold = 0, min.pct = 0.25,
                                     verbose = FALSE)
fwrite(cluster.markers_F3, "2.0_EpiHR_LNM/cluster.markers_F3.csv")

cluster.markers_F3 <- fread("2.0_EpiHR_LNM/cluster.markers_F3.csv",
                            data.table = FALSE)
DEG3 <- cluster.markers_F3[cluster.markers_F3$cluster %in% "pos", ]
DEG3$dif_pct <- DEG3$pct.1 - DEG3$pct.2
names(DEG3) <- c("pValue", "Log2FC", "pct.1", "pct.2",
                 "FDR", "cluster", "GeneName", "dif_pct")
rownames(DEG3) <- DEG3$GeneName
DEG3$Log2FC[DEG3$Log2FC < -2.5] <- -2.5
DEG3$FDR[DEG3$FDR < 1e-80] <- 1e-80

data <- add_regulate(DEG3, log2FC_name = "Log2FC",
                     fdr_name = "FDR", log2FC = 1, fdr = 0.05)
data$row <- rownames(data)

top15 <- cluster.markers_F3 %>% TOP_N(15, pct.1 = 0.2)
top15 <- top15[top15$avg_log2FC > 0 & top15$cluster %in% "pos", ]

p_volcano3 <- gradual_volcano1(
  data, x = "log2FoldChange", y = "padj",
  label = "row", label_number = 10,
  output = FALSE, custom_label = c(top15$gene, "CXCL1")
) + ggsci::scale_color_gsea() + ggsci::scale_fill_gsea()
ggsave("2.0_EpiHR_LNM/Volcano_F3_top15.pdf", p_volcano3,
       width = 6, height = 6)

# ---- 22. Mouse scRNA-seq analysis ------------------------------------------
scobj_mouse <- qs::qread("F3_mouse/scobj.final_mouse.qs")

# Figure 3A
p_mouse_anno <- CellDimPlot(srt = scobj_mouse, group.by = "anno",
                            show_stat = FALSE, pt.size = 1.5,
                            reduction = "umap", theme_use = "theme_blank")
ggsave("F3_mouse/DimPlot.anno.pdf", p_mouse_anno, width = 20, height = 10)

# Figure 3B
p_mouse_split <- CellDimPlot(srt = scobj_mouse, group.by = "anno",
                             show_stat = FALSE, split.by = "group",
                             pt.size = 1.5, reduction = "umap",
                             theme_use = "theme_blank")
ggsave("F3_mouse/DimPlot.group.anno.pdf", p_mouse_split,
       width = 20, height = 10)

in.dat <- scobj_mouse@meta.data
observe.data <- table(in.dat$group, in.dat$anno)
x <- chisq.test(observe.data)
R_oe <- x$observed / x$expected
plot_ro_e(R_oe, "F3_mouse/Roe.anno.pdf")

# Figure S3A,S3B: lineage markers
p_lineage <- FeatureDimPlot(
  scobj_mouse, show_stat = FALSE,
  features = c("Ptprc", "Epcam", "Col1a1", "Cd34",
               "Cd3e", "Klrb1c", "Cd19", "Cd14",
               "Cd68", "Fscn1", "Mcpt4"),
  xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 1.5, theme_use = "theme_blank"
)
ggsave("F3_mouse/FeatureDimPlot_lineage.pdf", p_lineage,
       width = 16, height = 12)

# ---- 23. T/NK cell analysis ------------------------------------------------
scobj_T <- qs::qread("F3_mouse/3.0_T_cell_main/scobj.anno.final.qs")

# Figure 3D
p_t_anno <- CellDimPlot(srt = scobj_T, group.by = "anno",
                        show_stat = FALSE, pt.size = 2,
                        reduction = "umap", theme_use = "theme_blank")
ggsave("F3_mouse/3.0_T_cell_main/DimPlot.anno.pdf", p_t_anno,
       width = 20, height = 10)

# Figure 3E: exhaustion markers
p_exhaust <- FeatureDimPlot(
  scobj_T, show_stat = FALSE,
  features = c("Tigit", "Pdcd1", "Havcr2", "Ctla4"),
  xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 3, theme_use = "theme_blank"
)
ggsave("F3_mouse/3.0_T_cell_main/FeatureDimPlot_Exhausted.marker.pdf",
       p_exhaust, width = 20, height = 20)

# Figure 3F: Ro/e for T cell subtypes
in.dat <- scobj_T@meta.data
observe.data <- table(in.dat$group, in.dat$anno)
x <- chisq.test(observe.data)
R_oe <- x$observed / x$expected
plot_ro_e(R_oe, "F3_mouse/3.0_T_cell_main/Roe.anno.pdf", transpose = TRUE)

# Figure 3G: markers and correlation with F3
Idents(scobj_T) <- "anno"
cluster.markers <- FindAllMarkers(scobj_T, assay = "RNA", only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25,
                                  verbose = FALSE)
write.csv(cluster.markers,
          "F3_mouse/3.0_T_cell_main/T.RNA_snn_res.0.2.markers.anoo.csv",
          row.names = FALSE)

cluster.markers1 <- cluster.markers
cluster.markers1$gene1 <- str_to_upper(cluster.markers1$gene)
cluster.markers1 <- cluster.markers1[
  cluster.markers1$pct.1 > 0.2 & cluster.markers1$pct.2 < 0.5, ]
top20 <- cluster.markers1 %>% TOP_N(20, pct.1 = 0.2)
epi_top20 <- lapply(split(top20, top20$cluster),
                    function(x) unique(x$gene1))

exp1 <- GSVA::gsva(rt_total, epi_top20, method = "ssgsea", parallel.sz = 8)
exp1_1 <- t(rbind(exp1, "F3" = rt_total["F3", ]))

corr <- round(cor(exp1_1), 2)
p.mat <- cor_pmat(exp1_1)
p_corr <- ggcorrplot(corr, type = "lower", outline.color = "black",
                     colors = c("#6D9EC1", "white", "#E46726"))
ggsave("F3_mouse/3.0_T_cell_main/corr.sum.pdf", p_corr,
       width = 6, height = 6)

p_scatter <- ggplot(data.frame(exp1_1), aes(F3, Treg)) +
  geom_point(col = "#984ea3") +
  geom_smooth(method = lm, se = TRUE, na.rm = TRUE,
              fullrange = TRUE, size = 2, col = "#fdc086") +
  geom_rug(col = "#7fc97f") +
  stat_cor(method = "pearson", digits = 3, size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave("F3_mouse/3.0_T_cell_main/F3_Treg_cor.pdf", p_scatter,
       width = 6, height = 6)

# Figure 3H: survival by Treg score
clin_total <- cbind(clin_total[, 1:19], exp1_1)
clin_total$F3 <- as.numeric(rt_total["F3", ])
clin_total$Treg_best <- ggsur_to_best(
  total = clin_total, gene = "Treg",
  OS_m = "OS_m", OS_status = "OS_status", minprop = 0.2
)[, 1]

pdf("F3_mouse/3.0_T_cell_main/Treg_best_survivals.pdf",width = 6,height = 6)
p_treg_surv <- ggsur_fun_cag(
  total = clin_total, gene = "Treg_best",
  color_palette1 = "d3", x_break = 36, conf.int = FALSE,
  OS_m = "OS_m", OS_status = "OS_status"
)
print(p_treg_surv, newpage = FALSE)
dev.off()



# Figure 4G: Foxp3 and Ccr6 FeaturePlot
p_foxp3 <- FeatureDimPlot(
  scobj_T, show_stat = FALSE,
  features = c("Foxp3", "Ccr6"),
  xlab = "UMAP1", ylab = "UMAP2",
  pt.size = 3, theme_use = "theme_blank"
)
ggsave("F3_mouse/3.0_T_cell_main/FeatureDimPlot_Foxp3_Ccr6.pdf",
       p_foxp3, width = 10, height = 10)

# Figure S3C: T cell UMAP split by group
p_t_split <- CellDimPlot(srt = scobj_T, group.by = "anno",
                         show_stat = FALSE, split.by = "group",
                         pt.size = 2, reduction = "umap",
                         theme_use = "theme_blank")
ggsave("F3_mouse/3.0_T_cell_main/DimPlot.group.anno.pdf", p_t_split,
       width = 20, height = 10)

# Figure S3D: DotPlot of selected markers
NK.cell.marker     <- c("Klrb1c", "Itga2", "Klrk1")
CD4.cell.marker    <- c("Cd4", "Foxp3", "Il2ra")
CD8.cell.marker    <- c("Cd8a", "Cd8b1")
Delta.cell.marker  <- c("Trdc", "Klrd1")
Exhausted.marker   <- c("Lag3", "Tigit", "Pdcd1", "Havcr2", "Ctla4", "Entpd1")
Cytotoxic.marker   <- c("Gzma", "Gzmb", "Ifng", "Nkg7")
markers <- unique(c(NK.cell.marker, CD4.cell.marker, CD8.cell.marker,
                    Delta.cell.marker, Exhausted.marker, Cytotoxic.marker))

Idents(scobj_T) <- "anno"
p_dot <- DotPlot(scobj_T, features = markers) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_gradientn(
    colours = c("#330066", "#336699", "#66CC66", "#FFCC33")
  ) +
  labs(x = NULL, y = NULL)
ggsave("F3_mouse/3.0_T_cell_main/DotPlot1.anno.pdf", p_dot,
       width = 6, height = 10)

# ---- 24. Session info ------------------------------------------------------
sessionInfo()