library(data.table)
library(dplyr)
library(ggplot2)
library(GSVA)
library(Biobase)
library(edgeR)
library(pheatmap)
library(tibble)
library(msigdbr)
library(cowplot)
library(tidyr)
library(scales)

set.seed(2025)

tcr_file <- "./7_.diversity.strict.resampled.txt"
rna_file <- "./aging_cohortTPM_nochi_forTIMER2.txt"
prot_file <- "./ZD_Lumos_ZD018_All_65Sample_Merge_Protein_UP1_Normalization.csv"
out_dir <- "./"

tcr <- fread(tcr_file, data.table = FALSE, check.names = FALSE)
rownames(tcr) <- tcr[, 1]
tcr <- tcr[, -1, drop = FALSE]

rna <- fread(rna_file, data.table = FALSE, check.names = FALSE)
rownames(rna) <- rna[, 1]
rna <- rna[, -1, drop = FALSE]

prot <- fread(prot_file, data.table = FALSE, check.names = FALSE)
prot <- prot[!duplicated(prot$Symbol), ]
rownames(prot) <- prot[, 2]
prot <- prot[, -c(1:3), drop = FALSE]

tcr_cols <- colnames(tcr)
pick_col <- function(pattern) {
  idx <- grep(pattern, tcr_cols, ignore.case = TRUE)
  if (length(idx) == 0) return(NULL)
  tcr_cols[idx]
}
tcr_metrics <- unique(c(pick_col("efron"), pick_col("chao1"), pick_col("shannon"), pick_col("diversity"), pick_col("inverseSimpson"), pick_col("observedDiversity")))
tcr_metrics <- tcr_metrics[!grepl("_std$", tcr_metrics)]
tcr_df <- tcr[, tcr_metrics, drop = FALSE]
tcr_df <- apply(tcr_df, 2, function(x) as.numeric(as.character(x)))
rownames(tcr_df) <- rownames(tcr)

common_samples <- Reduce(intersect, list(rownames(tcr_df), colnames(rna)))

tcr_df <- tcr_df[common_samples, , drop = FALSE]
rna <- rna[, common_samples, drop = FALSE]

rna_mat <- as.matrix(rna)
if (all(rna_mat == round(rna_mat)) && median(rna_mat, na.rm = TRUE) > 50) {
  dge <- DGEList(counts = rna_mat)
  cpm_mat <- cpm(dge, log = FALSE, prior.count = 1)
  expr <- log2(cpm_mat + 1)
} else {
  expr <- log2(rna_mat + 1)
}
rownames(expr) <- toupper(rownames(expr))

m_c7 <- msigdbr(species = "Homo sapiens", category = "C7")
m_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
gsets_c7 <- split(x = m_c7$gene_symbol, f = m_c7$gs_name)
gsets_go_bp <- split(x = m_c5$gene_symbol[m_c5$gs_subcat == "GO:BP"], 
                     f = m_c5$gs_name[m_c5$gs_subcat == "GO:BP"])

gsets_go_bp_sets <- gsets_go_bp[grepl("cytotox|cytotoxic|cytotoxicity|CD8|CTL", 
                                      names(gsets_go_bp), ignore.case = TRUE)]
gsets_go_bp_sets <- gsets_go_bp[grepl("exhausted", 
                                      names(gsets_go_bp), ignore.case = TRUE)]

selected_sets <- c("GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION",
  "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY"
  )
gsets_go_bp_sets <- gsets_go_bp[selected_sets]
gene_sets <- gsets_go_bp_sets

gset_upper <- lapply(gene_sets, toupper)
expr <- as.data.frame(expr)
expr <- expr[rowSums(expr) > 10, ]
expr_for_gsva <- expr
expr <- as.matrix(expr)
storage.mode(expr) <- "double"

gset_filtered <- gset_upper
if (length(gset_filtered) == 0) stop("No gene sets available.")

gsva_ver <- as.numeric(strsplit(as.character(packageVersion("GSVA")), "\\.")[[1]][1])

res_ssgsea <- NULL
try({
  ssg_par <- ssgseaParam(expr, gset_filtered,
                         minSize = 5,
                         maxSize = 500,
                         alpha = 0.25,
                         normalize = TRUE)
  res_ssgsea <- gsva(ssg_par, verbose = FALSE)
}, silent = TRUE)

if (is.null(res_ssgsea)) {
  res_ssgsea <- tryCatch({
    gsva(expr, gset_filtered, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE,
         min.sz = 5, max.sz = 500)
  }, error = function(e) {
    stop("GSVA failed: ", e$message)
  })
}

ssgsea_df <- t(res_ssgsea)

prot_mat <- as.matrix(prot)
if (median(prot_mat, na.rm = TRUE) > 50 && all(prot_mat == round(prot_mat))) {
  prot_expr <- log2(prot_mat + 1)
} else {
  prot_expr <- log2(prot_mat + 1)
}
rownames(prot_expr) <- toupper(rownames(prot_expr))

prot_targets <- c("GZMB", "PRF1", "GNLY", "GZMA", "NKG7", "CD8A", "CD3D")
prot_exist <- intersect(toupper(rownames(prot_expr)), prot_targets)

prot_df_sub <- if (length(prot_exist) > 0) t(prot_expr[prot_exist, , drop = FALSE]) else matrix(NA, nrow = ncol(prot_expr), ncol = 0)

prot_df_sub <- as.data.frame(prot_df_sub)
ssgsea_df <- as.data.frame(ssgsea_df)
ssgsea_df$ID <- rownames(ssgsea_df)
prot_df_sub$ID <- rownames(prot_df_sub)

predictor_df <- merge(ssgsea_df, prot_df_sub, by = 'ID', all.x = TRUE)

rownames(predictor_df) <- predictor_df$ID
predictor_df <- predictor_df[, -1]

results <- list()
for (metric in colnames(tcr_df)) {
  for (var in colnames(predictor_df)) {
    x <- as.numeric(tcr_df[, metric])
    y <- as.numeric(predictor_df[, var])
    ok <- which(!is.na(x) & !is.na(y))
    n_ok <- length(ok)
    if (n_ok >= 6) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      rho <- ct$estimate
      pval <- ct$p.value
    } else {
      rho <- NA
      pval <- NA
    }
    results[[length(results) + 1]] <- data.frame(
      metric = metric,
      variable = var,
      rho = as.numeric(rho),
      pval = as.numeric(pval),
      n = n_ok,
      stringsAsFactors = FALSE
    )
  }
}
res_df <- do.call(rbind, results)
res_df$padj <- p.adjust(res_df$pval, method = "BH")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)
annot_txt <- function(rho, p) sprintf("rho=%.3f\np=%.3g", rho, p)
res_df <- res_df[res_df$pval < 0.05, ]
res_df <- res_df[!is.na(res_df$metric), ]

for (i in seq_len(nrow(res_df))) {
  row <- res_df[i, ]
  if (is.na(row$rho)) next
  x <- as.numeric(tcr_df[, row$metric])
  y <- as.numeric(predictor_df[, row$variable])
  df_plt <- data.frame(sample = rownames(tcr_df), x = x, y = y)
  df_plt <- df_plt[!is.na(df_plt$x) & !is.na(df_plt$y), ]
  if (nrow(df_plt) < 6) next
  p <- ggplot(df_plt, aes(x = x, y = y)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = paste0(row$metric), y = paste0(row$variable),
         title = paste0("Spearman: ", row$metric, " vs ", row$variable),
         subtitle = annot_txt(row$rho, row$pval)) +
    theme_cowplot()
  ggsave(filename = file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", row$metric), "__vs__", gsub("[^A-Za-z0-9_\\-]", "_", row$variable), ".pdf")),
         plot = p, width = 6, height = 5)
}

res_df <- res_df[!is.na(res_df$rho), ]
df <- res_df

rho_mat_df <- df %>% dplyr::select(metric, variable, rho) %>%
  tidyr::pivot_wider(names_from = metric, values_from = rho)
rownames_mat <- rho_mat_df$variable
rho_mat <- as.matrix(rho_mat_df[, -1, drop = FALSE])
rownames(rho_mat) <- rownames_mat
rho_mat <- rho_mat[, -5]

get_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

sig_df <- df %>% dplyr::select(metric, variable, pval)
sig_wide <- tidyr::pivot_wider(sig_df, names_from = metric, values_from = pval)
sig_mat <- as.matrix(sig_wide[, -1, drop = FALSE])
rownames(sig_mat) <- sig_wide$variable
star_mat <- matrix("", nrow = nrow(sig_mat), ncol = ncol(sig_mat), dimnames = dimnames(sig_mat))
for (i in seq_len(nrow(sig_mat))) {
  for (j in seq_len(ncol(sig_mat))) {
    star_mat[i, j] <- get_star(as.numeric(sig_mat[i, j]))
  }
}

maxAbs <- max(abs(rho_mat), na.rm = TRUE)
breaks <- seq(-maxAbs, maxAbs, length.out = 101)
my_colors <- colorRampPalette(c("#039be5", "white", "lightcoral"))(length(breaks) - 1)

pdf('./TCR_metrics_vs_mRNA_protein_correlation_heatmap.pdf', width = 7.5, height = 5)
pheatmap(rho_mat,
         color = my_colors,
         breaks = breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = star_mat,
         number_color = "black",
         fontsize_number = 9,
         fontsize_row = 9, fontsize_col = 9,
         angle_col = 45,
         main = "Spearman rho: predictors vs TCR metrics")
dev.off()