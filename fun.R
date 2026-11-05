# =============================================================================
# fun.R -- Custom helper functions
# =============================================================================

# Required packages for these functions
# library(survival); library(survminer); library(ggplot2)
# library(ggrepel); library(latex2exp); library(dplyr); library(grid)

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- TOP_N: select top N markers per cluster --------------------------------
TOP_N <- function(x, n, pct.1 = 0.1, first = "avg_log2FC",
                  second = "p_val_adj", sig.padj = NULL,
                  fc.threshold = c(-0.25, 0.25)) {
  required_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2",
                     "p_val_adj", "gene", "cluster")
  if (!all(required_cols %in% names(x))) {
    message("be careful ! Required columns missing.")
    return(data.frame())
  }
  
  reordered_daf <- data.frame()
  
  if (!is.null(pct.1) && pct.1 > 0 && pct.1 < 1) {
    x <- x[x$pct.1 >= pct.1, ]
  }
  if (!is.null(sig.padj) && is.numeric(sig.padj)) {
    x <- x[x$p_val_adj <= sig.padj, ]
  }
  if (length(fc.threshold) == 2) {
    fc.threshold <- sort(fc.threshold)
    x <- x[(x$avg_log2FC < fc.threshold[1] | x$avg_log2FC > fc.threshold[2]), ]
  } else if (length(fc.threshold) == 1) {
    x <- x[x$avg_log2FC >= fc.threshold, ]
  }
  
  for (i in unique(x$cluster)) {
    sub_x <- x[x$cluster == i, ]
    if (n > nrow(sub_x)) {
      message("FBI warning, n < ", n)
      tmp <- dplyr::arrange(sub_x, desc(get(first)), desc(get(second)))
    } else {
      tmp <- dplyr::arrange(sub_x, desc(get(first)), desc(get(second)))[1:n, ]
    }
    reordered_daf <- rbind(reordered_daf, tmp)
  }
  return(reordered_daf)
}

# ---- cell_fun: custom heatmap cell drawing ---------------------------------
# NOTE: depends on global `plot.data` and `cols` defined before pheatmap call.
cell_fun <- function(j, i, x, y, width, height, fill) {
  if (plot.data[i, j] < as.numeric(names(cols)[1])) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = cols[1], fill = cols[1]))
    grid.text("\u00B1", x = x, y = y)
  } else if (plot.data[i, j] < as.numeric(names(cols)[2])) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = cols[2], fill = cols[2]))
    grid.text("+", x = x, y = y)
  } else if (plot.data[i, j] < as.numeric(names(cols)[3])) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = cols[3], fill = cols[3]))
    grid.text("++", x = x, y = y)
  } else {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = cols[4], fill = cols[4]))
    grid.text("+++", x = x, y = y)
  }
}

# ---- dif_test: t-test / Wilcoxon for a single gene --------------------------
dif_test <- function(clin, rt, group, gene) {
  rt_new <- data.frame(
    "group" = clin[, group],
    "gene"  = as.numeric(rt[gene, ])
  )
  t_dif    <- t.test(rt_new$gene ~ rt_new$group)
  wilcox_w <- wilcox.test(rt_new$gene ~ rt_new$group)
  
  data.frame(
    gene_name     = gene,
    z_test        = t_dif$statistic,
    pvalue_ttest  = t_dif$p.value,
    mean_high     = t_dif$estimate[2],
    mean_low      = t_dif$estimate[1],
    mean_dif      = t_dif$estimate[2] - t_dif$estimate[1],
    pvalue_wilcox = wilcox_w$p.value
  )
}

# ---- surv: Cox regression for a single gene --------------------------------
surv <- function(rt, clin, method = "continuous", time = "OS_m",
                 event = "OS_status", minprop = 0.1, gene) {
  rt1 <- data.frame(
    "time"  = clin[, time],
    "event" = clin[, event],
    "gene"  = as.numeric(rt[gene, ])
  )
  
  if (method == "continuous") {
    rt1$target_gene <- rt1$gene
  } else if (method == "median") {
    rt1$target_gene <- ifelse(rt1$gene > median(rt1$gene), 1, 0)
  } else if (method == "best") {
    suv_test <- surv_cutpoint(rt1, time = "time", event = "event",
                              variables = "gene", minprop = minprop)
    rt1$target_gene <- ifelse(rt1$gene > suv_test$cutpoint[1, 1], 1, 0)
  }
  
  cox        <- coxph(Surv(time, event) ~ target_gene, data = rt1)
  coxSummary <- summary(cox)
  
  data.frame(
    gene_name  = gene,
    z_survival = coxSummary$coefficients[, "z"],
    HR         = coxSummary$conf.int[, "exp(coef)"],
    HR.95L     = coxSummary$conf.int[, "lower .95"],
    HR.95H     = coxSummary$conf.int[, "upper .95"],
    pvalue_HR  = coxSummary$coefficients[, "Pr(>|z|)"]
  )
}

# ---- ggsur_to_best: dichotomize by optimal cutpoint -------------------------
ggsur_to_best <- function(total, gene, OS_m = "OS_m",
                          OS_status = "OS_status", minprop = 0.1) {
  res.cut <- surv_cutpoint(data = total, time = OS_m, minprop = minprop,
                           event = OS_status, variables = gene)
  a <- ifelse(total[, gene] > res.cut$cutpoint[1, 1], 1, 0)
  rt2 <- data.frame("a" = a)
  colnames(rt2) <- paste0(gene, "_best")
  return(rt2)
}

# ---- ggsur_fun_cag: Kaplan-Meier plot with HR annotation --------------------
ggsur_fun_cag <- function(total, gene, color_palette1, x_break,
                          conf.int = FALSE, legend.labs = c("low", "high"),
                          OS_m = "OS_m", OS_status = "OS_status") {
  rt_new <- data.frame(
    "time"   = total[, OS_m],
    "status" = total[, OS_status],
    "var"    = total[, gene]
  )
  fit1     <- survfit(Surv(time, status) ~ var, data = rt_new)
  fit_cox  <- coxph(Surv(time, status) ~ var, data = rt_new)
  coxSumm  <- summary(fit_cox)
  coxP     <- coxSumm$coefficients[, "Pr(>|z|)"]
  coxHR    <- format(round(coxSumm$conf.int[, "exp(coef)"], 2), nsmall = 2)
  coxCI    <- paste0(
    format(round(coxSumm$conf.int[, "lower .95"], 2), nsmall = 2), "-",
    format(round(coxSumm$conf.int[, "upper .95"], 2), nsmall = 2)
  )
  p_values <- paste0("HR: ", coxHR, ",", " 95% CI: ", coxCI,
                     "\n", "P = ", format(coxP, scientific = TRUE, digits = 3))
  
  ggsurv1 <- ggsurvplot(
    fit = fit1, data = rt_new,
    conf.int = conf.int,
    conf.int.alpha = 0.1,
    surv.median.line = "none",
    palette = color_palette1,
    xlim = c(0, round(max(fit1$time) / 6) * 6),
    ylim = c(0, 1.02),
    axes.offset = TRUE,
    xlab = "Time, mo",
    ylab = "Survival probability",
    break.time.by = x_break,
    surv.scale = "percent",
    censor.shape = "|", censor.size = 2,
    legend = c(0.8, 0.8),
    legend.labs = legend.labs,
    legend.title = gene,
    font.legend = c(16, "plain"),
    font.x = c(16, "plain", "black"),
    font.y = c(16, "plain", "black"),
    font.xtickslab = c(12, "plain", "black"),
    font.ytickslab = c(12, "plain", "black"),
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.fontsize = 4,
    risk.table.col = "strata",
    risk.table.y.text = TRUE,
    tables.y.text = FALSE,
    risk.table.height = 0.25,
    surv.plot.height = 0.75,
    title = "",
    font.main = c(16, "bold")
  )
  ggsurv1$table <- ggsurv1$table + labs(x = NULL, y = NULL)
  ggsurv1$plot <- ggsurv1$plot +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = 12, y = 0.25, label = p_values,
             cex = 5, vjust = 0.5, hjust = 0, fontface = 2, col = "black") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, round(max(fit1$time) / 6) * 6, x_break))
  return(ggsurv1)
}

# ---- gradual_volcano1: volcano plot with gradient colors -------------------
gradual_volcano1 <- function(data,
                             x = "log2FoldChange", y = "padj",
                             pointSizeRange = c(0.5, 4),
                             fills = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"),
                             colors = c("#17194e", "#68bfe7", "#f9ed36", "#a22f27", "#211f1f"),
                             x_lab = NULL, y_lab = NULL,
                             legend_title = NULL, legend_position = NULL,
                             add_line = TRUE, log2FC_cut = 1, FDR_cut = 0.05,
                             add_label = TRUE, label = "row",
                             label_number = 10, custom_label = NULL,
                             output = TRUE, filename = "volcano_plot") {
  colnames(data)[colnames(data) == x] <- "x"
  colnames(data)[colnames(data) == y] <- "y"
  colnames(data)[colnames(data) == label] <- "geneName"
  
  if (is.null(custom_label)) {
    data$label <- rep("", nrow(data))
    if (label_number != 0) {
      data$label[order(data$y)[1:label_number]] <-
        data$geneName[order(data$y)[1:label_number]]
    }
  } else {
    data$label <- rep("", nrow(data))
    data$label[match(custom_label, data$geneName)] <- custom_label
  }
  
  p <- ggplot(data, aes(x, -log10(y))) +
    geom_point(aes(size = -log10(y), fill = -log10(y), color = -log10(y)),
               shape = 21) +
    scale_fill_gradientn(values = seq(0, 1, 0.2), colors = fills) +
    scale_color_gradientn(values = seq(0, 1, 0.2), colors = colors) +
    scale_size_continuous(range = pointSizeRange) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.background = element_blank()) +
    guides(fill = guide_colourbar(title = legend_title %||% "-Log10_q-value"),
           color = "none", size = "none") +
    labs(x = x_lab %||% TeX("$Log_2 \\textit{FC}$"),
         y = y_lab %||% TeX("$-Log_{10} \\textit{FDR} $"))
  
  if (add_line) {
    p <- p +
      geom_vline(xintercept = c(-log2FC_cut, log2FC_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(FDR_cut), linetype = "dashed")
  }
  if (add_label) {
    p <- p +
      geom_text_repel(aes(label = label, color = -log10(y)),
                      size = 3, max.overlaps = 10000,
                      key_glyph = draw_key_point)
  }
  
  if (is.null(legend_position)) {
    p <- p + theme(legend.position = c(0.01, 0.7),
                   legend.justification = c(0, 1))
  } else if (legend_position == "UL") {
    p <- p + theme(legend.position = c(0.01, 0.99),
                   legend.justification = c(0, 1))
  } else if (legend_position == "UR") {
    p <- p + theme(legend.position = c(0.99, 0.99),
                   legend.justification = c(1, 1))
  } else if (legend_position == "DL") {
    p <- p + theme(legend.position = c(0.01, 0.01),
                   legend.justification = c(0, 0))
  } else {
    p <- p + theme(legend.position = c(0.99, 0.01),
                   legend.justification = c(1, 0))
  }
  
  if (output) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = 5, width = 6)
  }
  return(p)
}

# ---- getpansurdata: pan-cohort Cox results ---------------------------------
getpansurdata <- function(splitdata, method = "continuous", time = "OS_m",
                          event = "OS_status", minprop = 0.1, gene) {
  do.call(rbind, lapply(splitdata, function(rt) {
    rt1 <- data.frame(
      "time"  = rt[, time],
      "event" = rt[, event],
      "gene"  = as.numeric(rt[, gene])
    )
    if (method == "continuous") {
      rt1$target_gene <- rt1$gene
    } else if (method == "median") {
      rt1$target_gene <- ifelse(rt1$gene > median(rt1$gene), 1, 0)
    } else if (method == "best") {
      suv_test <- surv_cutpoint(rt1, time = "time", event = "event",
                                variables = "gene", minprop = minprop)
      rt1$target_gene <- ifelse(rt1$gene > suv_test$cutpoint[1, 1], 1, 0)
    }
    cox        <- coxph(Surv(time, event) ~ target_gene, data = rt1)
    coxSummary <- summary(cox)
    
    data.frame(
      gene_name  = gene,
      z_survival = coxSummary$coefficients[, "z"],
      HR         = coxSummary$conf.int[, "exp(coef)"],
      HR.95L     = coxSummary$conf.int[, "lower .95"],
      HR.95H     = coxSummary$conf.int[, "upper .95"],
      pvalue_HR  = coxSummary$coefficients[, "Pr(>|z|)"],
      type       = unique(rt$type)
    )
  }))
}