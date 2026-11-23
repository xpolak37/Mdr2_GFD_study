
suppressMessages(suppressWarnings({library(tidyverse)
  library(openxlsx)
  library(wesanderson)
  library(patchwork)
  library(reshape2)
  library(data.table)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(plotly)
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(UpSetR)
  library(gridExtra)
  library(ensembldb)
  library(EnsDb.Mmusculus.v79)
  library(ggpubr)
  library(vegan)
  library(emmeans)
  library(rstatix)
  library(kableExtra)
  library(knitr)
  library(DESeq2)
  library(ashr)
  library(tidyverse)
  library(data.table)
  library(openxlsx)
  library(vegan)
  library(emmeans)
  library(knitr)
  library(kableExtra)
  library(fgsea)
  library(msigdbr)
  library(patchwork)
  library(matrixStats)
  
}))


# Functions


# lm(expr ~ Group + Order), all genes at once
order_drift <- function(m, cd, formula = ~ Group + Order, term = "Order") {
  
  # build and check the design matrix
  X <- model.matrix(formula, cd)
  if (qr(X)$rank < ncol(X)) stop("Rank-deficient")
  n <- nrow(X); p <- ncol(X)
  
  # This computes (XᵀX)⁻¹, the core of the least-squares normal equations.
  XtXinv <- chol2inv(chol(crossprod(X)))
  dimnames(XtXinv) <- list(colnames(X), colnames(X))
  
  # solve all genes simultaneously and get test statistics
  beta   <- XtXinv %*% crossprod(X, t(m)) # the OLS estimate
  resid  <- t(m) - X %*% beta # observed minus fitted, per gene.
  sigma2 <- colSums(resid^2) / (n - p) # the residual variance for each gene, with the correct n − p degrees of freedom.
  
  # The actual test — is the Order slope ≠ 0?
  slope  <- beta[term, ]; se <- sqrt(XtXinv[term, term] * sigma2)
  tval   <- slope / se
  data.frame(gene_id = rownames(m), 
             slope = slope, 
             t = tval,
             p = 2 * pt(-abs(tval), n - p),
             padj = p.adjust(2 * pt(-abs(tval), n - p), "BH")) |>
    dplyr::left_join(sym, by = "gene_id")
}


summarise_all_res <- function(res, sym, name, txi, dds, 
                              sample_metadata,
                              groups, # e.g. c("SD_F","GFD_F")
                              group_vars = c("Diet", "Sex"), # build the label
                              experiment="3"
) {
  
  ## sample -> group label
  meta <- sample_metadata %>%
    unite("group", all_of(group_vars), sep = "_", remove = FALSE)
  
  sample_ids <- meta$SampleIAB    # the sample IDs 
  if (experiment=="4") sample_ids <- meta$SampleID_rna
  grp        <- meta$group
  names(grp) <- sample_ids
  
  ## helper: gene x sample matrix -> gene x group means, only for `groups`
  group_means <- function(mat, prefix) {
    mat <- mat[, sample_ids, drop = FALSE]  # align columns to metadata order
    out <- sapply(groups, function(lv)
      rowMeans(mat[, grp == lv, drop = FALSE]))
    colnames(out) <- paste0(prefix, "_", groups)
    as.data.frame(out) %>% rownames_to_column("gene_id")
  }
  
  ## results table 
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(sym, by = "gene_id") %>%
    arrange(padj)
  
  ## per-group means (only the requested groups) 
  tpm_grp  <- group_means(txi$abundance,                   "TPM")
  cnt_grp  <- group_means(counts(dds, normalized = FALSE), "count")
  norm_grp <- group_means(counts(dds, normalized = TRUE),  "normcount")
  
  ## join 
  res_df <- res_df %>%
    left_join(tpm_grp,  by = "gene_id") %>%
    left_join(cnt_grp,  by = "gene_id") %>%
    left_join(norm_grp, by = "gene_id")
  
  ## write 
  write.xlsx(res_df,paste0("../Results/",name,"_P",experiment,"_all.xlsx"),rownames=FALSE)
  return(res_df)
}


ggtheme <- function(p){
  p <- p + theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, 
                                      size = 1,linewidth = 1),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          axis.text = element_text(face = "bold",colour = "black"),
          axis.title = element_text(face = "bold",colour = "black"),
          panel.grid = element_blank(),
          legend.position = "right")
  return(p)
}




# volcano plot
# res        : a DESeqResults / lfcShrink object
# title      : panel title
# lfc_color  : LFC threshold
# clip       : x-axis clip; LFCs beyond ±clip are pulled to the edge
# xlab       : x-axis label (differs for the interaction panel)
# n_label    : how many top genes (by padj) to label
make_volcano <- function(res, title,
                         lfc_color = 1, clip = 10,
                         xlab = expression(log[2]~fold~change),
                         n_label = 12, alpha = 0.05) {
  df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    dplyr::left_join(sym, by = c("gene" = "gene_id")) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      lfc_plot = pmax(pmin(log2FoldChange, clip), -clip),   # clip extremes
      clipped  = abs(log2FoldChange) > clip,
      sig = case_when(
        padj < alpha & log2FoldChange >=  lfc_color ~ "Up",
        padj < alpha & log2FoldChange <= -lfc_color ~ "Down",
        TRUE                                        ~ "NS"
      )
    )
  
  top <- df %>% dplyr::filter(sig != "NS") %>% slice_min(padj, n = n_label)
  
  p1 <- ggplot(df, aes(lfc_plot, -log10(padj))) +
    geom_point(
      aes(color = sig, shape = clipped), 
      size = 1.6, alpha = 0.7) +
    scale_color_manual(
      values = c(Down = "#3B7DB3", NS = "grey80", Up = "#B33B3B"),
      name = "sig") +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17),   
      guide = "none") +
    geom_vline(xintercept = c(-lfc_color, lfc_color),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = top, aes(label = gene_name),
                             size = 3, max.overlaps = Inf,
                             min.segment.length = 0, seed = 42) +
    labs(title = title, 
         x = xlab, 
         y = expression(-log[10]~adjusted~italic(p))) +
    coord_cartesian(xlim = c(-clip, clip)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  
  return(p1)
}

make_volcano_publication <- function(res, title,
                         lfc_color = 1, clip = 10,
                         xlab = expression(log[2]~fold~change),
                         n_label = 12, alpha = 0.05) {
  df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    dplyr::left_join(sym, by = c("gene" = "gene_id")) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      lfc_plot = pmax(pmin(log2FoldChange, clip), -clip),   # clip extremes
      clipped  = abs(log2FoldChange) > clip,
      sig = case_when(
        padj < alpha & log2FoldChange >=  lfc_color ~ "Up",
        padj < alpha & log2FoldChange <= -lfc_color ~ "Down",
        TRUE                                        ~ "NS"
      )
    )
  
  n_sig <- sum(df$sig != "NS")
  count_label <- paste0(format(n_sig, big.mark = " "), " genes\n(females: ", 3, " genes)")
  # n_sig_f = number of significant genes in the female comparison; drop the 2nd line if not needed
  
  p2 <- ggplot(df, aes(lfc_plot, -log10(padj))) +
    # non-significant first, so coloured points sit on top
    geom_point(data = subset(df, sig == "NS"),
               colour = "grey85", size = 0.6) +
    geom_point(data = subset(df, sig != "NS"),
               aes(colour = sig, shape = clipped), size = 0.9) +
    scale_colour_manual(values = c(Down = "#2C4E80", Up = "#B03A2E")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    geom_vline(xintercept = c(-lfc_color, lfc_color),
               linetype = "dotted", colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(alpha),
               linetype = "dotted", colour = "grey40", linewidth = 0.3) +
    annotate("text", x = -clip, y = Inf, label = count_label,
             hjust = 0, vjust = 1.3, size = 2.6, lineheight = 0.9) +
    scale_x_continuous(breaks = seq(-clip, clip, by = 2)) +
    coord_cartesian(xlim = c(-clip, clip)) +
    labs(title = title, x = xlab,
         y = "-log10 adjusted p") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 10, face = "plain"),
      plot.tag        = element_text(size = 14, face = "bold"),
      axis.text       = element_text(colour = "black", size = 8),
      axis.line       = element_line(linewidth = 0.4),
      axis.ticks      = element_line(linewidth = 0.4)
    )
  
  
  return(p2)
}

pca_plot <- function(df,metadata,colors,perform_pca=TRUE,
                     axis=c("PC1","PC2")){
  axis_num <- c(as.numeric(regmatches(axis[1], regexpr ("\\d+" , axis[1] ))),
                as.numeric(regmatches(axis[2], regexpr ("\\d+" , axis[2] ))))
  
  if (perform_pca){
    log_df <- log2(df + 1)
    
    data.pca <- prcomp(log_df,scale. = TRUE,center = TRUE)
    # variance of each PC
    var_explained <- (data.pca$sdev)^2
    
    # proportion of variance explained
    prop_var <- var_explained / sum(var_explained) 
    
    # cumulative variance explained
    imp_vec <- cumsum(prop_var) *100
    
    x_lab = paste(axis[1], "(",round(prop_var[axis_num[1]],2) * 100,"%", ")", sep="")
    y_lab = paste(axis[2], "(",round(prop_var[axis_num[2]],2) * 100,"%", ")", sep="")
    pca_df  <- as.data.frame(data.pca$x)
  } else {
    x_lab <- axis[1]
    y_lab <- axis[2]
    pca_df <- df
  }
  
  pca_df <- merge(pca_df %>% rownames_to_column("SampleID"),
                  metadata,
                  by="SampleID")
  pca_df$SampleID <- gsub("S25_001_","S",pca_df$SampleID)
  p <- ggplot(pca_df) + 
    geom_point(aes(x=!!sym(axis[1]),y=!!sym(axis[2]),
                   color=Group,
                   text = paste0("Sample: ", SampleID,
                                 "<br>Group: ", Group)),size=4) +
    #geom_text_repel(aes(x=!!sym(axis[1]),y=!!sym(axis[2]),label=SampleID),size=3) + 
    xlab(x_lab) + 
    ylab(y_lab) + 
    scale_color_manual(values=colors)
  
  p <- ggtheme(p)
  return(p)
}

# tables
ktab <- function(df, caption, digits = 4) {
  df %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, digits))) %>%
    kable(caption = caption, booktabs = TRUE) %>%
    kable_styling(latex_options = "hold_position",
                  bootstrap_options = c("striped", "hover", "condensed"))
}

# DESeq2 results as data frame with symbols 
res_to_df <- function(res) {
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(sym, by = "gene_id")
}

# one-sided and equivalence p-values 
# expected = +1 (effect > 0) or -1 (effect < 0)
p_dir <- function(t, df, expected) {
  ifelse(expected > 0, pt(t, df, lower.tail = FALSE), pt(t, df, lower.tail = TRUE))
}

# is the effect inside [-margin, +margin]?
p_tost <- function(est, se, df, margin) {
  pmax(pt((est + margin) / se, df, lower.tail = FALSE),
       pt((est - margin) / se, df, lower.tail = TRUE))
}

# results log
# type: "dir"  = directional test, replicated if sign matches and p < alpha
#       "tost" = equivalence test, supported if p < alpha
#       "null" = assumption check, holds if p >= alpha
log_result <- function(id, # if of hypotesis (H1/H2..)
                       tier, # tier (primary..)
                       what, # name
                       test, # test name
                       type, # type of test (directional, tost, null)
                       expected = NA,
                       estimate = NA, 
                       p = NA, 
                       disc = NA, # discovery, statistics from
                       note = "") {
  results_log[[length(results_log) + 1]] <<- data.frame(
    id = id, tier = tier, what = what, test = test, type = type,
    expected = expected, P76_3 = disc, P76_4 = estimate, p = p, note = note
  )
  invisible(NULL)
}
