library(psych)
library(reshape2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(cowplot)
library(GSVA)
library(viper)
library(dorothea)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(motifmatchr)

write.table(clin_info, './clin_info.txt', col.names = T, row.names = F, sep = '\t', quote = F)

prot_norm <- read.csv('./ZD_Lumos_ZD018_All_65Sample_Merge_Protein_UP1_Normalization.csv', header = T, stringsAsFactors = F)
prot_norm_nochi <- prot_norm[, grep('XC', invert = T, colnames(prot_norm))]
prot_norm_nochi <- prot_norm_nochi[rowSums(!is.na(prot_norm_nochi[, 4:51])) > 0, ]
prot_norm_nochi <- prot_norm_nochi[!duplicated(prot_norm_nochi$Symbol), ]
rownames(prot_norm_nochi) <- prot_norm_nochi$Symbol
prot_norm_nochi <- prot_norm_nochi[, -c(1, 2, 3)]
all(clin_info$ID %in% colnames(prot_norm_nochi))
prot_norm_nochi <- prot_norm_nochi[, clin_info$ID]
prot_norm_nochi <- prot_norm_nochi[rowSums(!is.na(prot_norm_nochi)) >= 0.25 * 48, ]

exp_mat <- t(prot_norm_nochi)
age_df <- clin_info[, c(1, 2)]
rownames(age_df) <- age_df$ID
age_df <- age_df[rownames(exp_mat), ]
all(rownames(exp_mat) == rownames(age_df))
age_df <- age_df[2]
exp_mat <- as.matrix(exp_mat)
age_df <- as.matrix(age_df)
rownames(age_df) == rownames(exp_mat)
class(age_df[1, 1])
class(exp_mat[1, 1])

spearman_test <- corr.test(x = exp_mat, y = age_df, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df <- as.data.frame(spearman_test$r)
cor_df$id <- row.names(cor_df)
cor_df <- cor_df[, c(2, 1)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)
p_df <- as.data.frame(spearman_test$p)
p_df$id <- row.names(p_df)
p_df <- p_df[, c(2, 1)]
p_melt <- melt(p_df, value.name = 'V2')
p_melt <- unique(p_melt)
colnames(cor_melt) <- c('gene', 'age', 'cor')
cor_melt$p_adj <- p_melt$V2
res_spearman <- cor_melt
colnames(res_spearman)[1:2] <- c('gene', 'age')
res_filt <- res_spearman[res_spearman$p_adj < 0.05, ]
res_filt <- res_filt[!is.na(res_filt$cor), ]
age_prots <- res_filt
age_prots_up <- age_prots[age_prots$cor > 0, ]
age_prots_down <- age_prots[age_prots$cor < 0, ]

sf_heat <- prot_norm_for_heat[age_sf$gene[age_sf$cor < 0], ]
sf_heat <- sf_heat[, clin_info$ID]
sf_heat_scaled <- t(scale(t(sf_heat)))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ha <- HeatmapAnnotation(group = clin_info$zubie, col = list(group = c('A2' = '#9e9d24', 'A3' = '#673ab7', 'A4' = '#f44336')))
col_fun <- colorRamp2(c(-2, 0, 2), c("#00c853", "white", "#c51162"))
col_split <- c(rep(1, sum(clin_info$zubie == 'A2')), rep(2, sum(clin_info$zubie == 'A3')), rep(3, sum(clin_info$zubie == 'A4')))
pdf('./11SF_downregulation_with_age_heatmap.pdf', height = 7, width = 10)
Heatmap(sf_heat_scaled, col = col_fun, cluster_rows = T, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 0), row_names_side = "left", show_column_names = F, show_row_names = T, border = F, show_heatmap_legend = T, row_names_gp = gpar(fotsize = 58), heatmap_legend_param = list(title = "Absolute_score", at = c(-2, 0, 2), legend_height = unit(3, "cm"), legend_width = unit(0.4, "cm"), labels = c("Low", "median", "High"), title_position = "topcenter", direction = "vertical", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)), height = unit(5, "cm"), width = unit(20, "cm"), column_split = col_split, column_gap = unit(3, "mm"), top_annotation = ha)
dev.off()

age_prots_up_ord <- age_prots_up[order(age_prots_up$p_adj), ]
age_prots_down_ord <- age_prots_down[order(age_prots_down$p_adj), ]
prot_norm_for_heat <- prot_norm_nochi[c(age_prots_up_ord$gene, age_prots_down_ord$gene), ]
prot_norm_for_heat_scaled <- t(scale(t(prot_norm_for_heat)))
rowSums(prot_norm_for_heat_scaled, na.rm = T)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ha <- HeatmapAnnotation(group = clin_info$zubie, col = list(group = c('A2' = '#9e9d24', 'A3' = '#673ab7', 'A4' = '#f44336')))
col_fun <- colorRamp2(c(-1.3, 0, 1.3), c("#00c853", "white", "#c51162"))
col_split <- c(rep(1, sum(clin_info$zubie == 'A2')), rep(2, sum(clin_info$zubie == 'A3')), rep(3, sum(clin_info$zubie == 'A4')))
pdf('./349proteins_associated_with_age_heatmap.pdf', height = 7, width = 10)
Heatmap(prot_norm_for_heat_scaled, col = col_fun, cluster_rows = F, cluster_columns = FALSE, row_names_side = "left", show_column_names = F, show_row_names = T, border = F, show_heatmap_legend = T, row_names_gp = gpar(fotsize = 58), heatmap_legend_param = list(title = "Absolute_score", at = c(-2, 0, 2), legend_height = unit(3, "cm"), legend_width = unit(0.4, "cm"), labels = c("Low", "median", "High"), title_position = "topcenter", direction = "vertical", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)), height = unit(5, "cm"), width = unit(20, "cm"), column_split = col_split, column_gap = unit(3, "mm"), top_annotation = ha)
dev.off()

high_genes <- age_prots_up$gene
go_bp_up <- enrichGO(gene = high_genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_bp_up_res <- go_bp_up@result

low_genes <- age_prots_down$gene
go_bp_down <- enrichGO(gene = low_genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_bp_down_res <- go_bp_down@result

sel_desc <- c('generation of precursor metabolites and energy', 'electron transport chain', 'ATP metabolic process', 'mitochondrial transmembrane transport', 'monosaccharide metabolic process', 'cytoplasmic translation', 'mRNA splicing, via spliceosome', 'regulation of RNA splicing', 'RNA splicing')
go_res_comb <- rbind(go_bp_up_res, go_bp_down_res)
df_bar <- go_res_comb[go_res_comb$Description %in% sel_desc, ]
df_bar$minus_logp <- -log10(df_bar$pvalue)
df_bar$Description <- factor(df_bar$Description, levels = sel_desc)
theme_set(theme_cowplot())
pdf('./Proteins_associated_with_age_GO_BP_barplot.pdf', height = 8, width = 4)
ggplot(data = df_bar, mapping = aes(Description, minus_logp)) + geom_bar(position = "dodge", stat = "identity", width = 0.8) + theme(axis.text.x = element_text(angle = 90))
dev.off()

sel_desc2 <- c('mRNA splicing, via spliceosome', 'RNA splicing', 'T cell activation', 'T cell differentiation', 'cytoplasmic translation', 'RNA localization')
df_bar2 <- go_bp_down_res[go_bp_down_res$Description %in% sel_desc2, ]
df_bar2$minus_logp <- -log10(df_bar2$pvalue)
df_bar2$Description <- factor(df_bar2$Description, levels = sel_desc2)
theme_set(theme_cowplot())
pdf('./Proteins_downregulation_with_age_GO_BP_barplot.pdf', height = 8, width = 4)
ggplot(data = df_bar2, mapping = aes(Description, minus_logp)) + geom_bar(position = "dodge", stat = "identity", width = 0.8) + theme(axis.text.x = element_text(angle = 90))
dev.off()

age_genes_all <- read.table('./age_associated_allgene_noChi_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)
age_genes_all <- age_genes_all[age_genes_all$p_adjust < 0.05, ]
age_genes_all_up <- age_genes_all[age_genes_all$correlation > 0, ]
age_genes_all_down <- age_genes_all[age_genes_all$correlation < 0, ]
age_prots_up$gene[age_prots_up$gene %in% age_genes_all_up$gene]
age_prots_down$gene[age_prots_down$gene %in% age_genes_all_down$gene]

sel_markers <- c('B3GAT1', 'KLRG1', 'LAG3', 'TIGIT', 'CTLA4', 'HAVCR2')
clin_mrna <- read.table('./aging_cohortClini_nochi.txt', header = T, sep = '\t', stringsAsFactors = F)
rownames(clin_mrna) <- clin_mrna$SampleID
clin_mrna <- clin_mrna[3]
expr_mrna <- read.table('./aging_cohortTPM_nochi_forTIMER2.txt', header = T, sep = '\t', stringsAsFactors = F)
expr_mrna <- expr_mrna[sel_markers, ]
expr_mrna <- expr_mrna[, clin_mrna$SampleID]
expr_mrna <- as.data.frame(t(expr_mrna))
rownames(expr_mrna) == rownames(expr_mrna)
clin_mrna[1, 1]
expr_mrna[1, 1]

spearman_test2 <- corr.test(x = expr_mrna, y = clin_mrna, use = "pairwise", method = "spearman", adjust = "none", ci = F)
cor_df2 <- as.data.frame(spearman_test2$r)
cor_df2$id <- row.names(cor_df2)
cor_df2 <- cor_df2[, c(2, 1)]
cor_melt2 <- melt(cor_df2, value.name = 'V2')
cor_melt2 <- unique(cor_melt2)
p_df2 <- as.data.frame(spearman_test2$p)
p_df2$id <- row.names(p_df2)
p_df2 <- p_df2[, c(2, 1)]
p_melt2 <- melt(p_df2, value.name = 'V2')
p_melt2 <- unique(p_melt2)
colnames(cor_melt2) <- c('gene', 'age', 'cor')
cor_melt2$p_adj <- p_melt2$V2
res_spearman2 <- cor_melt2
colnames(res_spearman2)[1:2] <- c('gene', 'age')
df_plot_mrna <- cbind(clin_mrna, expr_mrna)
theme_set(theme_cowplot())
ggplot(df_plot_mrna, aes(x = age, y = B3GAT1)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")

imm_gene_list <- readRDS('./ImmuneCellGeneList_ZhangQi.rds')
write.table(clin_info, './clin_info.txt', col.names = T, row.names = F, sep = '\t', quote = F)
prot_for_gsva <- prot_norm_nochi
prot_for_gsva <- prot_for_gsva[!duplicated(prot_for_gsva$Symbol), ]
rownames(prot_for_gsva) <- prot_for_gsva$Symbol
prot_for_gsva <- prot_for_gsva[, clin_info$ID]
prot_for_gsva <- log2(prot_for_gsva + 1)
prot_for_gsva <- as.data.frame(prot_for_gsva)
prot_for_gsva <- as.matrix(prot_for_gsva)
dim(prot_for_gsva)
prot_for_gsva[1, 1]
gsva_res <- gsva(prot_for_gsva, imm_gene_list, verbose = TRUE, method = "ssgsea")
rownames(gsva_res) == clin_info$ID
gsva_res <- as.data.frame(t(gsva_res))
rownames(gsva_res) == clin_info$ID
age_df2 <- clin_info
rownames(age_df2) <- age_df2$ID
age_df2 <- age_df2[2]
spearman_test3 <- corr.test(x = gsva_res, y = age_df2, use = "pairwise", method = "spearman", adjust = "none", ci = F)
cor_df3 <- as.data.frame(spearman_test3$r)
cor_df3$id <- row.names(cor_df3)
cor_df3 <- cor_df3[, c(2, 1)]
cor_melt3 <- melt(cor_df3, value.name = 'V2')
cor_melt3 <- unique(cor_melt3)
p_df3 <- as.data.frame(spearman_test3$p)
p_df3$id <- row.names(p_df3)
p_df3 <- p_df3[, c(2, 1)]
p_melt3 <- melt(p_df3, value.name = 'V2')
p_melt3 <- unique(p_melt3)
colnames(cor_melt3) <- c('gene', 'age', 'cor')
cor_melt3$p_adj <- p_melt3$V2
res_spearman3 <- cor_melt3
colnames(res_spearman3)[1:2] <- c('gene', 'age')

gsva_age <- gsva_res
gsva_age$age <- clin_info$age
theme_set(theme_cowplot())
pdf('./Cytotoxic_T_GSVA_VS_Age.pdf', height = 4, width = 4)
ggplot(gsva_age, aes(x = age, y = Cytotoxic_T)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()
pdf('./Tcm_GSVA_VS_Age.pdf', height = 4, width = 4)
ggplot(gsva_age, aes(x = age, y = Tcm)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()

gsva_age_group <- gsva_age
gsva_age_group$group <- clin_info$zubie
gene_box <- function(gene = "CCL14", group = "group", comparisons = my_comps, data) {
  my_theme <- theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid = element_blank())
  p <- ggboxplot(data, x = group, y = gene, ylab = gene, xlab = group, fill = group, palette = c('#999927', '#594096', '#e64537'), ggtheme = theme_set(my_theme))
  p + stat_compare_means(comparisons = my_comps)
}
name <- colnames(gsva_age_group)[1:14]
my_comps <- list(c('A2', 'A3'), c('A2', 'A4'), c('A3', 'A4'))
pdf('./massSpectrometry_GSVA_VS_AGE_box_plot.pdf', height = 10, width = 15)
p <- lapply(X = name, gene_box, data = gsva_age_group)
library(gridExtra)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()

sen_list <- readRDS('./SenMayo_list_forGSVA.rds')
cell_age_up <- readRDS('./CellAge_induceAging_list_forGSVA.rds')
cell_age_down <- readRDS('./CellAge_inhibitAging_list_forGSVA.rds')
imm_port <- readRDS('./immPort_17_immunePathway_list_forGSVA.rds')
aging_immune_list <- c(sen_list, cell_age_up, cell_age_down, imm_port)

gsva_res2 <- gsva(prot_for_gsva, aging_immune_list, verbose = TRUE, method = "ssgsea")
gsva_res2 <- as.data.frame(t(gsva_res2))
rownames(gsva_res2) == clin_info$ID
age_df3 <- clin_info
rownames(age_df3) <- age_df3$ID
age_df3 <- age_df3[2]
all(rownames(gsva_res2) == rownames(age_df3))

spearman_test4 <- corr.test(x = gsva_res2, y = age_df3, use = "pairwise", method = "spearman", adjust = "none", ci = F)
cor_df4 <- as.data.frame(spearman_test4$r)
cor_df4$id <- row.names(cor_df4)
cor_df4 <- cor_df4[, c(2, 1)]
cor_melt4 <- melt(cor_df4, value.name = 'V2')
cor_melt4 <- unique(cor_melt4)
p_df4 <- as.data.frame(spearman_test4$p)
p_df4$id <- row.names(p_df4)
p_df4 <- p_df4[, c(2, 1)]
p_melt4 <- melt(p_df4, value.name = 'V2')
p_melt4 <- unique(p_melt4)
colnames(cor_melt4) <- c('gene', 'age', 'cor')
cor_melt4$p_adj <- p_melt4$V2
res_spearman4 <- cor_melt4
colnames(res_spearman4)[1:2] <- c('gene', 'age')

gsva_age2 <- gsva_res2
gsva_age2$age <- clin_info$age
theme_set(theme_cowplot())
pdf('./CellAge_induceAging_GSVA_VS_Age.pdf', height = 4, width = 4)
ggplot(gsva_age2, aes(x = age, y = CellAge_induceAging)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()
pdf('./TCRsignalingPathway_GSVA_VS_Age.pdf', height = 4, width = 4)
ggplot(gsva_age2, aes(x = age, y = TCRsignalingPathway)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()

age_lncrna <- read.table('./age_associated_lncRNA_noChi_LNCipedia_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)
expr_lncrna <- read.table('./aging_cohortTPM_noChi_filtered_numeric_lncRNA.txt', header = T, sep = '\t', stringsAsFactors = F)
age_lncrna_expr <- expr_lncrna[age_lncrna$gene, ]
age_prots_expr <- prot_norm_nochi[age_prots$gene, ]
inter_samples <- colnames(age_prots_expr)[colnames(age_prots_expr) %in% colnames(age_lncrna_expr)]
age_lncrna_expr <- age_lncrna_expr[, inter_samples]
age_prots_expr <- age_prots_expr[, inter_samples]
age_lncrna_expr <- t(age_lncrna_expr)
age_prots_expr <- t(age_prots_expr)
all(rownames(age_lncrna_expr) == rownames(age_prots_expr))

spearman_test5 <- corr.test(x = age_lncrna_expr, y = age_prots_expr, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df5 <- as.data.frame(spearman_test5$r)
cor_df5$id <- row.names(cor_df5)
colnames(cor_df5)
ncol(cor_df5)
cor_df5 <- cor_df5[, c(350, 1:349)]
cor_melt5 <- melt(cor_df5, value.name = 'V2')
cor_melt5 <- unique(cor_melt5)
p_df5 <- as.data.frame(spearman_test5$p)
p_df5$id <- row.names(p_df5)
p_df5 <- p_df5[, c(350, 1:349)]
p_melt5 <- melt(p_df5, value.name = 'V2')
p_melt5 <- unique(p_melt5)
colnames(cor_melt5) <- c('gene', 'age', 'cor')
cor_melt5$p_adj <- p_melt5$V2
res_spearman5 <- cor_melt5
colnames(res_spearman5)[1:2] <- c('gene', 'age')
res_filt5 <- res_spearman5[res_spearman5$p_adj < 0.05, ]
res_filt5 <- res_filt5[!is.na(res_filt5$cor), ]
res_filt5_08 <- res_filt5[abs(res_filt5$cor) > 0.8, ]
length(unique(res_filt5_08$gene))
length(unique(res_filt5_08$age))

age_prots$prot_group <- 'up_with_age'
age_prots$prot_group[age_prots$cor < 0] <- 'down_with_age'
age_lncrna$lnc_group <- 'up_with_age'
age_lncrna$lnc_group[age_lncrna$correlation < 0] <- 'down_with_age'
age_prots_merge <- age_prots[, c(1, 5)]
age_lncrna_merge <- age_lncrna[, c(1, 7)]
res_filt5_08 <- merge(res_filt5_08, age_prots_merge, by.x = 'age', by.y = 'gene', all.x = T)
res_filt5_08 <- merge(res_filt5_08, age_lncrna_merge, by.x = 'gene', by.y = 'gene', all.x = T)
res_filt5_08$cor_group <- 'positive'
res_filt5_08$cor_group[res_filt5_08$cor < 0] <- 'negative'
table(res_filt5_08$lnc_group, res_filt5_08$cor_group)

res_filt5 <- merge(res_filt5, age_prots_merge, by.x = 'age', by.y = 'gene', all.x = T)
res_filt5 <- merge(res_filt5, age_lncrna_merge, by.x = 'gene', by.y = 'gene', all.x = T)
res_filt5$cor_group <- 'positive'
res_filt5$cor_group[res_filt5$cor < 0] <- 'negative'
table(res_filt5$lnc_group, res_filt5$cor_group)

res_filt5_07 <- res_filt5[abs(res_filt5$cor) > 0.7, ]
table(res_filt5_07$lnc_group, res_filt5_07$cor_group)

down_prots <- res_filt5_07$age[res_filt5_07$prot_group %in% 'down_with_age']
down_prots <- unique(as.character(down_prots))
go_bp_down_prots <- enrichGO(gene = down_prots, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_bp_down_prots_res <- go_bp_down_prots@result

down_sen <- down_prots[down_prots %in% sen_list$SenMayo]
down_up_aging <- down_prots[down_prots %in% cell_age_up$CellAge_induceAging]
down_down_aging <- down_prots[down_prots %in% cell_age_down$CellAge_inhibitAging]

expr_rna <- read.table('./aging_cohortTPM_nochi_forTIMER2.txt', header = T, sep = '\t', stringsAsFactors = F)
expr_rna <- expr_rna[(rowSums(expr_rna) > 10), ]
age_genes <- read.table('./age_associated_ProteinCodingGene_noChi_Gencode_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)
age_genes <- age_genes$gene

regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
tf_targets <- regulons %>% split(.$tf) %>% map(~ unique(.$target))
bg_genes <- rownames(expr_rna) %>% unique()
age_genes_int <- intersect(age_genes, bg_genes)
res_enrich <- tibble(tf = names(tf_targets)) %>%
  rowwise() %>%
  mutate(targets = list(intersect(tf_targets[[tf]], bg_genes)),
         k = sum(targets %in% age_genes_int),
         M = length(bg_genes),
         K = length(unique(tf_targets[[tf]])),
         n = length(age_genes_int),
         p = phyper(k - 1, K, M - K, n, lower.tail = FALSE),
         odds = (k / (n - k + 1e-6)) / ((K - k + 1e-6) / (M - K - (n - k) + 1e-6))) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  arrange(padj)
top_tfs <- res_enrich %>% filter(padj < 0.05) %>% arrange(padj)
plot_df <- res_enrich %>% arrange(padj)
plot_df <- plot_df[1:20, ]
p_enrich <- ggplot(plot_df, aes(x = reorder(tf, -log10(padj + 1e-300)), y = -log10(padj + 1e-300))) +
  geom_point(aes(size = k, color = odds)) +
  coord_flip() +
  labs(x = "", y = "-log10(BH adj P)", title = "TF regulon enrichment in age-associated genes") +
  scale_size_continuous(name = "Overlap count") +
  scale_color_viridis_c(name = "Odds ratio") +
  theme_bw(base_size = 14)
pdf("./TF_enrichment_dotplot.pdf", width = 5.5, height = 4.5)
print(p_enrich)
dev.off()

regulons_df <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
norm_mor <- function(x) {
  x <- tolower(as.character(x))
  x[x %in% c("1", "+1", "+", "activation", "activator", "act")] <- "1"
  x[x %in% c("-1", "-", "repression", "repressor", "rep")] <- "-1"
  x[!(x %in% c("1", "-1"))] <- "1"
  as.numeric(x)
}
regulons_df <- regulons_df %>% mutate(tfmode = norm_mor(mor))
regulons_df <- regulons_df %>% mutate(likelihood = 1.0)
regulons_split <- split(regulons_df, regulons_df$tf)
regulons_list <- lapply(regulons_split, function(sub) {
  targets <- as.character(sub$target)
  tfmode_vec <- as.numeric(sub$tfmode)
  names(tfmode_vec) <- targets
  lik_vec <- as.numeric(sub$likelihood)
  names(lik_vec) <- targets
  if (any(duplicated(names(tfmode_vec)))) {
    uniq <- unique(names(tfmode_vec))
    tfmode_vec2 <- setNames(numeric(length(uniq)), uniq)
    lik_vec2 <- setNames(numeric(length(uniq)), uniq)
    for (g in uniq) {
      idx <- which(names(tfmode_vec) == g)
      svals <- tfmode_vec[idx]
      svals <- svals[!is.na(svals)]
      if (length(svals) == 0) ssign <- 1 else ssign <- ifelse(mean(svals) >= 0, 1, -1)
      tfmode_vec2[g] <- ssign
      lik_vec2[g] <- mean(lik_vec[idx], na.rm = TRUE)
    }
    tfmode_vec <- tfmode_vec2
    lik_vec <- lik_vec2
  }
  list(tfmode = tfmode_vec, likelihood = lik_vec)
})
names(regulons_list) <- names(regulons_split)
expr_genes <- rownames(expr_rna)
reg_targets <- unique(unlist(lapply(regulons_list, function(x) names(x$tfmode))))
common_genes <- intersect(expr_genes, reg_targets)
regulons_list_filt <- lapply(regulons_list, function(z) {
  keep <- intersect(names(z$tfmode), expr_genes)
  if (length(keep) == 0) return(NULL)
  list(tfmode = z$tfmode[keep], likelihood = z$likelihood[keep])
})
regulons_list_filt <- regulons_list_filt[!sapply(regulons_list_filt, is.null)]
set.seed(2025)
tf_activity <- viper::viper(eset = as.matrix(expr_rna), regulon = regulons_list_filt, verbose = TRUE)
dim(tf_activity)
age_vec <- read.table('./aging_cohortClini_nochi.txt', header = T, sep = '\t', stringsAsFactors = F)
rownames(age_vec) <- age_vec$SampleID
age_vec <- age_vec[colnames(tf_activity), ]
all(age_vec$SampleID == colnames(tf_activity))
age_vec <- as.vector(age_vec$age)
tf_corr <- apply(tf_activity, 1, function(x) {
  ct <- cor.test(x, age_vec, method = "spearman", exact = FALSE)
  c(rho = as.numeric(ct$estimate), p = ct$p.value)
})
tf_corr_df <- data.frame(tf = rownames(tf_activity),
                         rho = as.numeric(tf_corr["rho", ]),
                         p = as.numeric(tf_corr["p", ]),
                         stringsAsFactors = FALSE)
tf_corr_df$padj <- p.adjust(tf_corr_df$p, method = "BH")
tf_corr_df <- tf_corr_df[order(tf_corr_df$padj), ]

top_tfs2 <- top_tfs %>% head(4)
top_tfs2 <- top_tfs2$tf
top_tfs2 <- top_tfs2[top_tfs2 %in% rownames(tf_activity)]
act_plot_df <- as.data.frame(t(tf_activity[top_tfs2, , drop = FALSE]))
act_plot_df$age <- age_vec
act_plot_df_long <- pivot_longer(act_plot_df, cols = -age, names_to = "tf", values_to = "activity")
ann_list <- lapply(top_tfs2, function(tf) {
  act <- as.numeric(tf_activity[tf, ])
  ct <- suppressWarnings(cor.test(act, age_vec, method = "spearman", exact = FALSE))
  rho <- as.numeric(ct$estimate)
  pval <- ct$p.value
  p_str <- ifelse(is.na(pval), "NA", ifelse(pval < 0.001, "<0.001", sprintf("%.3g", pval)))
  lab <- paste0("rho = ", sprintf("%.2f", rho), "\nP = ", p_str)
  age_min <- min(age_vec, na.rm = TRUE)
  age_max <- max(age_vec, na.rm = TRUE)
  x_pos <- age_min + 0.78 * (age_max - age_min)
  y_min <- min(act, na.rm = TRUE)
  y_max <- max(act, na.rm = TRUE)
  if (is.na(y_min) || is.na(y_max) || y_max == y_min) {
    y_pos <- y_max
  } else {
    y_pos <- y_min + 0.88 * (y_max - y_min)
  }
  data.frame(tf = tf, rho = rho, p = pval, label = lab, x = x_pos, y = y_pos, stringsAsFactors = FALSE)
})
ann_df <- do.call(rbind, ann_list)
act_plot_df_long$tf <- factor(act_plot_df_long$tf, levels = top_tfs2)
ann_df$tf <- factor(ann_df$tf, levels = top_tfs2)
p_activity <- ggplot(act_plot_df_long, aes(x = age, y = activity)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE, size = 0.6) +
  facet_wrap(~tf, scales = "free_y", ncol = 2) +
  geom_text(data = ann_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.2, fontface = "bold") +
  labs(x = "Age", y = "Inferred TF activity (VIPER)", title = "Top TF activities vs Age") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5))
cairo_pdf("./Top_TF_activity_vs_Age_前4大调控最多年龄相关mRNA的TF(有p和spearman).pdf", width = 6, height = 5)
print(p_activity)
dev.off()

sel_tfs_heat <- top_tfs %>% head(21)
sel_tfs_heat <- sel_tfs_heat$tf
sel_tfs_heat <- sel_tfs_heat[sel_tfs_heat %in% rownames(tf_activity)]
mat <- tf_activity[sel_tfs_heat, , drop = FALSE]
ord <- order(age_vec)
annotation_col <- HeatmapAnnotation(Age = age_vec[ord])
ht <- Heatmap(mat[, ord], name = "TF activity",
              show_row_names = TRUE, show_column_names = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE,
              top_annotation = annotation_col,
              col = colorRamp2(c(min(mat), 0, max(mat)), c("blue", "white", "red")))
cairo_pdf("./TF_activity_heatmap.pdf", width = 5, height = 4)
draw(ht, heatmap_legend_side = "right")
dev.off()

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
tss_df <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "transcription_start_site", "strand"),
                filters = "hgnc_symbol", values = bg_genes, mart = mart)
keep_chroms <- c(as.character(1:22), "X", "Y", "MT", "M")
tss_df <- tss_df %>% filter(chromosome_name %in% keep_chroms)
tss_df <- tss_df %>%
  mutate(chrom = case_when(
    chromosome_name %in% c("MT", "M") ~ "chrM",
    TRUE ~ paste0("chr", chromosome_name)
  ))
tss_df <- tss_df %>% filter(!is.na(transcription_start_site))
tss_df <- tss_df %>% group_by(hgnc_symbol, chrom) %>%
  dplyr::slice(1) %>% ungroup()
upstream <- 2000
downstream <- 200
tss_df <- tss_df %>%
  mutate(strand2 = ifelse(as.numeric(strand) == 1, "+", "-"),
         tss = as.integer(transcription_start_site))
gr <- GRanges(
  seqnames = tss_df$chrom,
  ranges = IRanges(
    start = ifelse(tss_df$strand2 == "+", tss_df$tss - upstream, tss_df$tss - downstream),
    end = ifelse(tss_df$strand2 == "+", tss_df$tss + downstream, tss_df$tss + upstream)
  ),
  strand = tss_df$strand2,
  gene = tss_df$hgnc_symbol
)
start(gr) <- pmax(start(gr), 1L)
bs_seqlvls <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
canonical_bs <- paste0("chr", c(1:22, "X", "Y"))
canonical_bs <- c(canonical_bs, "chrM")
keep <- intersect(seqlevels(gr), canonical_bs)
gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
seqlen <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
end(gr) <- pmin(end(gr), seqlen[as.character(seqnames(gr))])
valid_idx <- which(start(gr) <= end(gr))
if (length(valid_idx) < length(gr)) {
  warning("去除 ", length(gr) - length(valid_idx), " 个无效 promoter 区间（start > end）。")
  gr <- gr[valid_idx]
}
prom_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
names(prom_seqs) <- mcols(gr)$gene

pfm_list <- TFBSTools::getMatrixSet(JASPAR2022, opts = list(species = 9606, matrixtype = "PWM"))
matches <- matchMotifs(pfm_list, prom_seqs, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores")
motif_counts <- motifCounts(matches)

safe_try_A <- function() {
  mm <- tryCatch({
    matchMotifs(pfm_list, subject = prom_seqs, genome = BSgenome.Hsapiens.UCSC.hg38)
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(mm)) return(NULL)
  mm_mat <- tryCatch({
    motifmatchr::motifMatches(mm)
  }, error = function(e) {
    return(NULL)
  })
  if (!is.null(mm_mat)) {
    mm_mat <- as(mm_mat, "matrix")
    motif_presence <- mm_mat
    motif_counts_binary <- ifelse(motif_presence, 1L, 0L)
    rownames(motif_counts_binary) <- rownames(mm_mat)
    colnames(motif_counts_binary) <- colnames(mm_mat)
    return(list(method = "A", motif_counts = motif_counts_binary, motif_presence = motif_presence))
  } else {
    return(NULL)
  }
}
resA <- safe_try_A()
if (!is.null(resA)) {
  motif_presence_mat <- resA$motif_presence
  motif_counts_mat <- resA$motif_counts
}
rownames(motif_presence_mat) <- names(prom_seqs)
promoter_genes <- names(prom_seqs)
age_idx <- promoter_genes %in% age_genes_int

if (exists("motif_presence_mat")) {
  mat_raw <- motif_presence_mat
} else if (exists("motif_counts_mat")) {
  mat_raw <- (motif_counts_mat > 0)
} else stop("找不到 motif_presence_mat 或 motif_counts_mat，请先运行 motif scanning 步骤。")
mat_raw <- as.matrix(mat_raw)
row_hits <- sum(rownames(mat_raw) %in% age_genes_int, na.rm = TRUE)
col_hits <- sum(colnames(mat_raw) %in% age_genes_int, na.rm = TRUE)
if (col_hits > row_hits) {
  motif_by_prom <- mat_raw
  motif_ids <- rownames(motif_by_prom)
  promoter_names <- colnames(motif_by_prom)
} else if (row_hits > col_hits) {
  motif_by_prom <- t(mat_raw)
  motif_ids <- rownames(motif_by_prom)
  promoter_names <- colnames(motif_by_prom)
} else {
  row_is_motif_like <- all(grepl("^MA\\d+", rownames(mat_raw)))
  col_is_motif_like <- all(grepl("^MA\\d+", colnames(mat_raw)))
  if (col_is_motif_like & !row_is_motif_like) {
    motif_by_prom <- t(mat_raw)
    motif_ids <- rownames(motif_by_prom)
    promoter_names <- colnames(motif_by_prom)
  } else if (row_is_motif_like & !col_is_motif_like) {
    motif_by_prom <- mat_raw
    motif_ids <- rownames(motif_by_prom)
    promoter_names <- colnames(motif_by_prom)
  } else {
    stop("无法自动判定矩阵方向：请确认矩阵的行/列名或手动设置 promoter_names 与 motif_ids。")
  }
}

M <- length(promoter_names)
age_idx <- promoter_names %in% age_genes_int
n_age <- sum(age_idx)
res_list <- lapply(seq_len(nrow(motif_by_prom)), function(i) {
  pres <- as.logical(motif_by_prom[i, ])
  a <- sum(pres & age_idx, na.rm = TRUE)
  b <- sum(pres & !age_idx, na.rm = TRUE)
  c_val <- n_age - a
  d_val <- (M - n_age) - b
  ft <- tryCatch(fisher.test(matrix(c(a, b, c_val, d_val), nrow = 2), alternative = "greater"),
                 error = function(e) list(estimate = NA, p.value = NA))
  or <- if (!is.null(ft$estimate)) as.numeric(ft$estimate) else NA
  data.frame(motif = motif_ids[i],
             n_present_age = a,
             n_present_bg = b,
             total_with_motif = sum(pres, na.rm = TRUE),
             odds = or,
             p = ft$p.value,
             stringsAsFactors = FALSE)
})
motif_enrich_df <- bind_rows(res_list)
motif_enrich_df$padj <- p.adjust(motif_enrich_df$p, method = "BH")
motif_enrich_df <- motif_enrich_df %>% arrange(padj)
if (exists("pfm_list")) {
  map_motif_to_tf <- function(id) {
    if (id %in% names(pfm_list)) {
      pfm <- pfm_list[[id]]
      name <- tryCatch(slot(pfm, "name"), error = function(e) NA_character_)
      if (is.na(name) || nchar(name) == 0) name <- tryCatch(slot(pfm, "ID"), error = function(e) NA_character_)
      if (is.na(name) || nchar(name) == 0) {
        tags <- tryCatch(slot(pfm, "tags"), error = function(e) NULL)
        if (!is.null(tags) && is.list(tags) && "geneSymbol" %in% names(tags)) name <- tags$geneSymbol
      }
      return(ifelse(is.null(name) || is.na(name), NA_character_, as.character(name)))
    } else return(NA_character_)
  }
  motif_enrich_df$TF_name <- sapply(motif_enrich_df$motif, map_motif_to_tf)
} else motif_enrich_df$TF_name <- NA_character_

topn <- 20
plot_df <- motif_enrich_df %>% dplyr::slice(1:topn) %>%
  mutate(neglog10p = -log10(p + 1e-320),
         label = ifelse(!is.na(TF_name), paste0(motif, " (", TF_name, ")"), motif))
p <- ggplot(plot_df, aes(x = reorder(label, neglog10p), y = neglog10p, size = n_present_age, color = odds)) +
  geom_point() + coord_flip() +
  labs(x = "", y = "-log10(P)", title = "Top motif enrichment in promoters of age-associated genes") +
  theme_bw(base_size = 12) +
  scale_size_continuous(name = "Count in age genes") +
  scale_color_continuous(name = "Odds ratio")
cairo_pdf("./motif_enrichment_dotplot_top20.pdf", width = 8, height = 6)
print(p)
dev.off()

out_pdf <- "./MYC_motif_enrichment_plot.pdf"
pdf_w <- 5
pdf_h <- 4
show_padj_label <- F
if (!exists("motif_enrich_df")) stop("找不到 motif_enrich_df，请先确保该表在当前环境中。")
df <- motif_enrich_df
req_cols <- c("motif")
if (!all(c("n_present_age", "n_present_bg") %in% colnames(df))) {
  stop("motif_enrich_df 必须包含 'n_present_age' 和 'n_present_bg' 两列（各 motif 在 age 基因和非-age 背景中的计数）。")
}
if (!("p" %in% colnames(df)) && !("padj" %in% colnames(df))) {
  stop("motif_enrich_df 必须包含 'p' 或 'padj' 列。")
}
myc_rows <- c()
if ("TF_name" %in% colnames(df)) {
  myc_rows <- which(grepl("\\bMYC\\b|MYCN|c-Myc|myc|MAX", df$TF_name, ignore.case = TRUE))
}
if (length(myc_rows) == 0) {
  myc_rows <- which(grepl("MYC|MYCN|c-Myc|myc|MAX", df$motif, ignore.case = TRUE))
}
if (length(myc_rows) == 0) {
  stop("在 motif_enrich_df 中未自动检测到与 MYC 相关的 motif。若有具体 motif IDs（例如 'MAxxxx.x'），请把它们放到向量 myc_ids 并重试。")
}
df_myc <- df[myc_rows, , drop = FALSE]
rownames(df_myc) <- NULL

n_age <- NA
M <- NA
if (all(c("n", "M") %in% colnames(df))) {
  n_age <- df$n[1]
  M <- df$M[1]
}
if (is.na(n_age) || is.na(M)) {
  if (exists("motif_by_prom") && is.matrix(motif_by_prom)) {
    promoter_names <- colnames(motif_by_prom)
    n_age <- length(intersect(promoter_names, age_genes_int))
    M <- length(promoter_names)
  } else if (exists("promoter_names") && is.vector(promoter_names)) {
    n_age <- length(intersect(promoter_names, age_genes_int))
    M <- length(promoter_names)
  } else if (exists("age_genes_int") && exists("prom_seqs") && !is.null(names(prom_seqs))) {
    promoter_names <- names(prom_seqs)
    n_age <- length(intersect(promoter_names, age_genes_int))
    M <- length(promoter_names)
  }
}
res_list <- lapply(seq_len(nrow(df_myc)), function(i) {
  row <- df_myc[i, ]
  a <- as.numeric(row$n_present_age)
  b <- as.numeric(row$n_present_bg)
  K <- a + b
  if (!is.na(n_age) && !is.na(M)) {
    c_val <- n_age - a
    d_val <- (M - n_age) - b
    tbl <- matrix(c(a, b, c_val, d_val), nrow = 2)
    ft <- tryCatch(fisher.test(tbl, alternative = "greater"), error = function(e) NULL)
    if (!is.null(ft) && !is.na(ft$estimate)) {
      or_est <- as.numeric(ft$estimate)
      ci_low <- as.numeric(ft$conf.int[1])
      ci_high <- as.numeric(ft$conf.int[2])
      pval <- as.numeric(ft$p.value)
    } else {
      a2 <- a + 0.5
      b2 <- b + 0.5
      c2 <- c_val + 0.5
      d2 <- d_val + 0.5
      or_est <- (a2 * d2) / (b2 * c2)
      se_log <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
      ci_low <- exp(log(or_est) - 1.96 * se_log)
      ci_high <- exp(log(or_est) + 1.96 * se_log)
      pval <- if ("p" %in% colnames(row)) as.numeric(row$p) else NA_real_
    }
  } else {
    a2 <- a + 0.5
    b2 <- b + 0.5
    or_est <- (a2) / (b2)
    ci_low <- NA_real_
    ci_high <- NA_real_
    pval <- if ("p" %in% colnames(row)) as.numeric(row$p) else NA_real_
  }
  padj_val <- if ("padj" %in% colnames(row)) as.numeric(row$padj) else NA_real_
  data.frame(motif = as.character(row$motif),
             TF_name = if ("TF_name" %in% colnames(row)) as.character(row$TF_name) else NA_character_,
             n_present_age = a,
             n_present_bg = b,
             OR = or_est,
             CI_low = ci_low,
             CI_high = ci_high,
             p = pval,
             padj = padj_val,
             stringsAsFactors = FALSE)
})
res_df <- bind_rows(res_list)
if (all(is.na(res_df$padj)) && any(!is.na(res_df$p))) {
  res_df$padj <- p.adjust(res_df$p, method = "BH")
}
res_plot <- res_df %>% filter(!is.na(OR) & OR > 0)
if (nrow(res_plot) == 0) {
  warning("无可绘制的 MYC motif OR（OR 全为 NA 或 <=0）。已写出表格但未生成图。")
  write.csv(res_df, "./MYC_motif_enrichment_table.csv", row.names = FALSE)
} else {
  res_plot <- res_plot %>% mutate(label = ifelse(!is.na(TF_name) & nzchar(TF_name),
                                                 paste0(motif, " (", TF_name, ")"), motif),
                                  logOR = log10(OR),
                                  neglog10padj = -log10(padj + 1e-300))
  p <- ggplot(res_plot, aes(x = reorder(label, OR), y = OR)) +
    geom_point(aes(size = n_present_age, color = neglog10padj), stroke = 0.6) +
    coord_flip() +
    scale_y_log10(labels = scales::label_number(scale = 1, accuracy = 0.01)) +
    scale_color_viridis_c(name = "-log10(padj)") +
    scale_size_continuous(name = "Count in age genes", range = c(2, 6)) +
    labs(x = "", y = "Odds Ratio (log scale)",
         title = "Enrichment of MYC-related motif(s) in promoters of age-associated genes") +
    theme_minimal()
  if (show_padj_label) {
    p <- p + geom_text(aes(label = ifelse(!is.na(padj) & padj < 0.1, sprintf("padj=%.3g", padj), "")),
                       hjust = -0.05, size = 3, color = "black", data = res_plot, na.rm = TRUE)
  }
  cairo_pdf(out_pdf, width = pdf_w, height = pdf_h)
  print(p)
  dev.off()
  write.csv(res_df, "./MYC_motif_enrichment_table.csv", row.names = FALSE)
}