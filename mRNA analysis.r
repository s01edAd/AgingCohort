library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(psych)
library(reshape2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Biobase)
library(DESeq2)
library(GSVA)
library(ggsci)

aging_tpm <- read.table('./aging_cohortTPM.txt', header = T, sep = '\t', stringsAsFactors = F)
aging_clini <- read.csv('./RNAgroup_aging_from_ZhangQi.csv', header = T, stringsAsFactors = F)
tpm_nochi <- aging_tpm[, -c(1:8)]
clini_nochi <- aging_clini[-c(1:8), ]
rownames(clini_nochi) <- 1:nrow(clini_nochi)
write.table(clini_nochi, './aging_cohortClini_nochi.txt', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(tpm_nochi, './aging_cohortTPM_nochi_forTIMER2.txt', col.names = T, row.names = T, quote = F, sep = '\t')

timer2_res <- read.csv('./TIMER2_output_estimation_matrix.csv', header = T, stringsAsFactors = F)
timer2_res <- timer2_res[grep('.*CIBERSORT$', timer2_res$cell_type), ]

cibersort_res <- read.csv('./CIBERSORTx_RelativeMode_output.csv', header = T, stringsAsFactors = F)
rownames(cibersort_res) <- cibersort_res$Mixture
cibersort_res <- cibersort_res[, -1]
cibersort_res <- as.data.frame(t(cibersort_res))
cibersort_res <- cibersort_res[1:22, ]
rowSums(cibersort_res)
cibersort_res <- cibersort_res[, clini_nochi$SampleID]

colnames(cibersort_res) == clini_nochi$SampleID
cibersort_kw <- cibersort_res
cibersort_kw <- cibersort_kw[rowSums(cibersort_kw) != 0, ]

kw_res <- matrix(nrow = nrow(cibersort_kw), ncol = 1)
kw_res <- as.data.frame(kw_res)
for (i in 1:nrow(cibersort_kw)) {
  kw_test <- kruskal.test(as.numeric(cibersort_kw[i, ]), clini_nochi$group)
  p <- kw_test$p.value
  kw_res[i, 1] <- p
}
colnames(kw_res) <- 'KW_P_value'
rownames(kw_res) <- rownames(cibersort_kw)
kw_res$Immune_cell <- rownames(kw_res)
write.table(kw_res, './kruskal_test_for_CIBERSORTxResult_between3group.txt', col.names = T, row.names = F, sep = '\t', quote = F)

cibersort_scaled <- t(scale(t(cibersort_res)))
cibersort_scaled <- as.matrix(cibersort_scaled)
rowSums(cibersort_scaled)
cibersort_scaled[is.nan(cibersort_scaled)] <- 0
clini_nochi$SampleID == colnames(cibersort_scaled)
clini_nochi$SampleID %in% colnames(cibersort_scaled)
colnames(cibersort_scaled) %in% clini_nochi$SampleID
cibersort_scaled <- cibersort_scaled[, clini_nochi$SampleID]
clini_nochi$SampleID == colnames(cibersort_scaled)

col_fun <- colorRamp2(c(-2, 0, 2), c("#00c853", "white", "#c51162"))
column_split <- c(rep(1, 16), rep(2, 22), rep(3, 29))
df <- data.frame(type = c(rep(1, 16), rep(2, 22), rep(3, 29)))
col_anno <- columnAnnotation(df = df, col = list(type = c('1' = '#9e9d24', '2' = '#673ab7', '3' = '#f44336')))
pdf('./CIBERSORTx_RelativeMode_AgingCohort_nochi.pdf', height = 12, width = 12)
Heatmap(cibersort_scaled, col = col_fun,
        cluster_rows = T, cluster_columns = F,
        rect_gp = gpar(col = "gray", lwd = 0.1), row_names_side = "left",
        show_column_names = F, show_row_names = T, border = FALSE, show_heatmap_legend = TRUE,
        row_names_gp = gpar(fotsize = 58),
        heatmap_legend_param = list(title = "Absolute_score", at = c(-2, 0, 2),
                                    legend_height = unit(3, "cm"), legend_width = unit(0.4, "cm"), labels = c("Low", "median", "High"),
                                    title_position = "topcenter", direction = "vertical",
                                    labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
        height = unit(10, "cm"), width = unit(10, "cm"),
        column_split = column_split,
        column_gap = unit(3, "mm"),
        top_annotation = col_anno)
dev.off()

my_comps <- list(c('Youth', 'Middle age'),
                 c('Middle age', 'The elderly'),
                 c('Youth', 'The elderly'))

gene_box <- function(gene = "CCL14", group = "Group", comparisons = my_comps, data) {
  my_theme <- theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.grid = element_blank())
  p <- ggboxplot(data, x = group, y = gene,
                 ylab = gene,
                 xlab = group,
                 fill = group,
                 palette = c('#999927', '#594096', '#e64537'),
                 ggtheme = theme_set(my_theme))
  p + stat_compare_means(comparisons = my_comps, method = 'wilcox.test')
}

cibersort_box <- cibersort_res
cibersort_box <- as.data.frame(t(cibersort_box))
rownames(cibersort_box) == clini_nochi$SampleID
cibersort_box$group <- clini_nochi$group

name <- colnames(cibersort_box)[1:22]
pdf(file = "./CIBERSORTx_result_batch_boxplot.pdf", width = 15, height = 20)
p <- lapply(X = name, gene_box, group = 'group', data = cibersort_box)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()

write.table(clini_nochi, './aging_cohortClini_nochi.txt', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(tpm_nochi, './aging_cohortTPM_nochi_forTIMER2.txt', col.names = T, row.names = T, quote = F, sep = '\t')
rownames(clini_nochi) <- clini_nochi$SampleID

tpm_filter <- tpm_nochi[(rowSums(tpm_nochi) >= ncol(tpm_nochi)), ]
exp_mat <- t(tpm_filter)
age_mat <- clini_nochi[3]
exp_mat <- exp_mat[rownames(age_mat), ]
rownames(exp_mat)
all(rownames(exp_mat) == rownames(age_mat))
exp_mat <- as.matrix(exp_mat)
age_mat <- as.matrix(age_mat)
spearman_res <- corr.test(x = exp_mat, y = age_mat, use = "pairwise",
                          method = "spearman",
                          adjust = "BH", ci = F)

cor_df <- as.data.frame(spearman_res$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(2, 1)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)

p_adj <- as.data.frame(spearman_res$p)
p_adj$ASE_id <- row.names(p_adj)
p_adj <- p_adj[, c(2, 1)]
p_adj_melt <- melt(p_adj, value.name = 'V2')
p_adj_melt <- unique(p_adj_melt)

colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
cor_melt$p_adjust <- p_adj_melt$V2
exp_age_spearman <- cor_melt
colnames(exp_age_spearman)[1:2] <- c('gene', 'age')
write.table(exp_age_spearman, './age_associated_allgene_noChi_SpearmanResult.txt', col.names = T, row.names = F, sep = '\t', quote = F)

exp_age_spearman_005 <- exp_age_spearman[exp_age_spearman$p_adjust < 0.05, ]
age_assoc_allgene <- exp_age_spearman_005

age_assoc_allgene$ENSEMBL <- mapIds(org.Hs.eg.db, keys = age_assoc_allgene$gene, column = 'ENSEMBL', keytype = 'SYMBOL')
sum(is.na(age_assoc_allgene$ENSEMBL))
age_assoc_allgene$ENSEMBL_SYMBOL <- age_assoc_allgene$ENSEMBL
age_assoc_allgene$ENSEMBL_SYMBOL[is.na(age_assoc_allgene$ENSEMBL_SYMBOL)] <- age_assoc_allgene$gene[is.na(age_assoc_allgene$ENSEMBL_SYMBOL)]

protein_coding_list <- read.delim('./ProteinCodingGene_list_from_gencodeV41.txt', header = T, sep = '\t', stringsAsFactors = F)

age_assoc_protein <- age_assoc_allgene[age_assoc_allgene$gene %in% protein_coding_list$gene_name, ]
write.table(age_assoc_protein, './age_associated_ProteinCodingGene_noChi_Gencode_SpearmanResult.txt', col.names = T, row.names = F, sep = '\t')

write.table(exp_age_spearman, './age_associated_allgene_noChi_SpearmanResult.txt', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(age_assoc_protein, './age_associated_ProteinCodingGene_noChi_Gencode_SpearmanResult.txt', col.names = T, row.names = F, sep = '\t')

age_assoc_allgene <- read.table('./age_associated_allgene_noChi_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)
age_assoc_allgene <- age_assoc_allgene[age_assoc_allgene$p_adjust < 0.05, ]

tmp1 <- as.data.frame((age_mat))
tmp2 <- as.data.frame((exp_mat))
tmp2$age <- tmp1$age

ggplot(tmp2, aes(x = age, y = ABAT)) + geom_point() + geom_smooth(method = "lm", color = "black", fill = "lightgray")

age_assoc_high <- age_assoc_allgene[age_assoc_allgene$correlation > 0, ]
high_exp <- age_assoc_high$gene

GO_en_BP <- enrichGO(gene = high_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './Genes_upregulation_with_age_GO_BP.csv', row.names = F, col.names = T, quote = F)

age_assoc_low <- age_assoc_allgene[age_assoc_allgene$correlation < 0, ]
low_exp <- age_assoc_low$gene

GO_en_BP <- enrichGO(gene = low_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './Genes_downregulation_with_age_GO_BP.csv', row.names = F, col.names = T, quote = F)

clini_nochi_ordered <- clini_nochi[order(clini_nochi$age), ]
tpm_nochi_ordered <- tpm_nochi[, clini_nochi_ordered$SampleID]
tpm_nochi_ordered_log2 <- log2(tpm_nochi_ordered + 1)

age_list <- unique(clini_nochi_ordered$age)
tpm_avg <- matrix(nrow = nrow(tpm_nochi_ordered_log2))
tpm_avg <- as.data.frame(tpm_avg)

for (i in age_list) {
  selected_sample <- clini_nochi_ordered$SampleID[clini_nochi_ordered$age == i]
  selected_sample <- as.character(selected_sample)
  if (length(selected_sample) == 1) {
    onecol <- as.data.frame(tpm_nochi_ordered_log2[, selected_sample])
    colnames(onecol) <- paste0('age_', i)
    df <- onecol
  } else {
    morecol <- as.data.frame(tpm_nochi_ordered_log2[, selected_sample])
    avg_col <- as.data.frame(apply(morecol, MARGIN = 1, FUN = mean))
    colnames(avg_col) <- paste0('age_', i)
    df <- avg_col
  }
  tpm_avg <- cbind(tpm_avg, df)
}
tpm_avg <- tpm_avg[, -1]
tpm_avg <- as.matrix(tpm_avg)

avg_pheno <- data.frame(Sample = colnames(tpm_avg))
avg_pheno$Age <- gsub('age_', '', avg_pheno$Sample)
avg_pheno$index <- 0:(nrow(avg_pheno) - 1)
avg_pheno <- avg_pheno[, c(3, 1, 2)]
colnames(avg_pheno) <- c('index', 'label', 'time')
avg_pheno$time <- as.numeric(avg_pheno$time)
rownames(avg_pheno) <- avg_pheno$label

tpm_avg <- read.table('./AgingCohort_TPM_noChi_log2_Average.txt', header = T, sep = '\t', stringsAsFactors = F)

selected_col <- c('age_23', 'age_30', 'age_40', 'age_50', 'age_60', 'age_70')
cm <- clusterData(exp = tpm_avg[, selected_col], cluster.method = "mfuzz", cluster.num = 8)

cm_filter <- cm
cm_filter$long.res <- cm_filter$long.res[cm_filter$long.res$membership >= 0.7, ]
cm_filter$wide.res <- cm_filter$wide.res[cm_filter$wide.res$membership >= 0.7, ]

table(cm_filter$long.res$cluster)

visCluster(object = cm_filter, plot.type = "line")
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"), add.mline = FALSE)

clini_nochi_ordered <- clini_nochi[order(clini_nochi$age), ]
tpm_nochi_ordered <- tpm_nochi[, clini_nochi_ordered$SampleID]

All_point <- list(which(clini_nochi_ordered$age >= 23 & clini_nochi_ordered$age <= 25),
                  which(clini_nochi_ordered$age >= 27 & clini_nochi_ordered$age <= 30),
                  which(clini_nochi_ordered$age >= 36 & clini_nochi_ordered$age <= 47),
                  which(clini_nochi_ordered$age >= 49 & clini_nochi_ordered$age <= 54),
                  which(clini_nochi_ordered$age >= 60 & clini_nochi_ordered$age <= 65),
                  which(clini_nochi_ordered$age >= 66 & clini_nochi_ordered$age <= 71))

tpm_selected <- matrix(nrow = nrow(tpm_nochi_ordered), ncol = length(All_point))
tpm_selected <- as.data.frame(tpm_selected)
rownames(tpm_selected) <- rownames(tpm_nochi_ordered)

for (i in 1:length(All_point)) {
  selected_point <- All_point[[i]]
  if (length(selected_point) > 1) {
    point_df <- rowMeans(tpm_nochi_ordered[, selected_point])
    tpm_selected[, i] <- point_df
  } else {
    point_df <- tpm_nochi_ordered[, selected_point]
    tpm_selected[, i] <- point_df
  }
}
colnames(tpm_selected) <- paste0('point', 1:length(All_point))

tpm_selected <- tpm_selected[rowMeans(tpm_selected) > 1, ]
tpm_selected <- log2(tpm_selected + 1)

cm <- clusterData(exp = tpm_selected, cluster.method = "mfuzz", cluster.num = 15)
visCluster(object = cm, plot.type = "line")

cm_filter <- cm
cm_filter$long.res <- cm_filter$long.res[cm_filter$long.res$membership >= 0.7, ]
cm_filter$wide.res <- cm_filter$wide.res[cm_filter$wide.res$membership >= 0.7, ]
table(cm_filter$long.res$cluster)
visCluster(object = cm_filter, plot.type = "line")
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"), add.mline = FALSE)

All_point <- list(which(clini_nochi_ordered$age >= 23 & clini_nochi_ordered$age <= 25),
                  which(clini_nochi_ordered$age >= 27 & clini_nochi_ordered$age <= 30),
                  which(clini_nochi_ordered$age >= 47 & clini_nochi_ordered$age <= 50),
                  which(clini_nochi_ordered$age >= 51 & clini_nochi_ordered$age <= 54),
                  which(clini_nochi_ordered$age >= 60 & clini_nochi_ordered$age <= 65),
                  which(clini_nochi_ordered$age >= 66 & clini_nochi_ordered$age <= 71))

tpm_selected <- matrix(nrow = nrow(tpm_nochi_ordered), ncol = length(All_point))
tpm_selected <- as.data.frame(tpm_selected)
rownames(tpm_selected) <- rownames(tpm_nochi_ordered)

for (i in 1:length(All_point)) {
  selected_point <- All_point[[i]]
  if (length(selected_point) > 1) {
    point_df <- rowMeans(tpm_nochi_ordered[, selected_point])
    tpm_selected[, i] <- point_df
  } else {
    point_df <- tpm_nochi_ordered[, selected_point]
    tpm_selected[, i] <- point_df
  }
}
colnames(tpm_selected) <- paste0('point', 1:length(All_point))

tpm_selected <- tpm_selected[rowMeans(tpm_selected) > 1, ]
tpm_selected <- log2(tpm_selected + 1)

cm <- clusterData(exp = tpm_selected, cluster.method = "mfuzz", cluster.num = 10)
visCluster(object = cm, plot.type = "line", add.mline = FALSE)

cm_filter <- cm
cm_filter$long.res <- cm_filter$long.res[cm_filter$long.res$membership >= 0.7, ]
cm_filter$wide.res <- cm_filter$wide.res[cm_filter$wide.res$membership >= 0.7, ]
table(cm_filter$long.res$cluster)
visCluster(object = cm_filter, plot.type = "line")
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(object = cm, plot.type = "line", ms.col = c("green", "orange", "red"), add.mline = FALSE)

genes_high <- age_assoc_high$gene
genes_low <- age_assoc_low$gene

genes_high_mat <- tpm_avg[genes_high, ]
genes_high_mat <- genes_high_mat[rowSums(is.na(genes_high_mat)) == 0, ]

genes_low_mat <- tpm_avg[genes_low, ]
genes_low_mat <- genes_low_mat[rowSums(is.na(genes_low_mat)) == 0, ]

cm <- clusterData(exp = genes_high_mat, cluster.method = "mfuzz", cluster.num = 9)
visCluster(object = cm, plot.type = "line", add.mline = FALSE)

All_point <- list(which(clini_nochi_ordered$age > 20 & clini_nochi_ordered$age <= 40),
                  which(clini_nochi_ordered$age > 40 & clini_nochi_ordered$age <= 60),
                  which(clini_nochi_ordered$age > 60 & clini_nochi_ordered$age <= 80))

tpm_selected <- matrix(nrow = nrow(tpm_nochi_ordered), ncol = length(All_point))
tpm_selected <- as.data.frame(tpm_selected)
rownames(tpm_selected) <- rownames(tpm_nochi_ordered)

for (i in 1:length(All_point)) {
  selected_point <- All_point[[i]]
  if (length(selected_point) > 1) {
    point_df <- rowMeans(tpm_nochi_ordered[, selected_point])
    tpm_selected[, i] <- point_df
  } else {
    point_df <- tpm_nochi_ordered[, selected_point]
    tpm_selected[, i] <- point_df
  }
}
colnames(tpm_selected) <- paste0('point', 1:length(All_point))

tpm_selected <- tpm_selected[rowMeans(tpm_selected) > 1, ]
tpm_selected <- log2(tpm_selected + 1)

tpm_high <- tpm_selected[genes_high, ]
tpm_high <- tpm_high[rowSums(is.na(tpm_high)) == 0, ]

tpm_low <- tpm_selected[genes_low, ]
tpm_low <- tpm_low[rowSums(is.na(tpm_low)) == 0, ]

cm <- clusterData(exp = tpm_high, cluster.method = "mfuzz", cluster.num = 4)
visCluster(object = cm, plot.type = "line", add.mline = FALSE)
table(cm$wide.res$cluster)

cm_filter <- cm
cm_filter$long.res <- cm_filter$long.res[cm_filter$long.res$membership >= 0.4, ]
cm_filter$wide.res <- cm_filter$wide.res[cm_filter$wide.res$membership >= 0.4, ]
table(cm_filter$wide.res$cluster)

pdf('./age_positive_genes_time_cluster_line.pdf', height = 4, width = 16)
visCluster(object = cm_filter, plot.type = "line", add.mline = F)
dev.off()

GO_res_high <- data.frame()
for (i in unique(cm_filter$wide.res$cluster)) {
  genes <- cm_filter$wide.res$gene[cm_filter$wide.res$cluster == i]
  GO_en_BP <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  GO_en_BP_res <- GO_en_BP@result
  GO_en_BP_res$cluster <- i
  GO_res_high <- rbind(GO_res_high, GO_en_BP_res)
}

GOterm_selected <- data.frame(id = c(rep('C1', 4), rep('C2', 2), rep('C3', 3), rep('C4', 2)),
                              term = c('cellular response to external stimulus',
                                       'positive regulation of interleukin-1 production',
                                       'defense response to virus',
                                       'T cell activation involved in immune response',
                                       'positive regulation of I-kappaB kinase/NF-kappaB signaling',
                                       'positive regulation of defense response',
                                       'positive regulation of I-kappaB kinase/NF-kappaB signaling',
                                       'myeloid cell development',
                                       'myeloid cell differentiation',
                                       'negative regulation of alpha-beta T cell activation',
                                       'negative regulation of alpha-beta T cell differentiation'))

pdf('./age_positive_genes_time_cluster_heatmap.pdf', height = 8, width = 15)
visCluster(object = cm_filter, plot.type = "both", ht.col = c("#00c853", "white", "#c51162"),
           add.box = T, add.line = T, line.side = "left", show_row_dend = F,
           annoTerm.data = GOterm_selected, boxcol = ggsci::pal_npg()(3))
dev.off()

cm <- clusterData(exp = tpm_low, cluster.method = "mfuzz", cluster.num = 4)
visCluster(object = cm, plot.type = "line", add.mline = FALSE)
table(cm$wide.res$cluster)

cm_filter <- cm
cm_filter$long.res <- cm_filter$long.res[cm_filter$long.res$membership >= 0.5, ]
cm_filter$wide.res <- cm_filter$wide.res[cm_filter$wide.res$membership >= 0.5, ]
table(cm_filter$wide.res$cluster)

pdf('./age_negative_genes_time_cluster_line.pdf', height = 4, width = 16)
visCluster(object = cm_filter, plot.type = "line", add.mline = F)
dev.off()

GO_res_low <- data.frame()
for (i in unique(cm_filter$wide.res$cluster)) {
  genes <- cm_filter$wide.res$gene[cm_filter$wide.res$cluster == i]
  GO_en_BP <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  GO_en_BP_res <- GO_en_BP@result
  GO_en_BP_res$cluster <- i
  GO_res_low <- rbind(GO_res_low, GO_en_BP_res)
}

GOterm_selected <- data.frame(id = c(rep('C1', 1), rep('C2', 3), rep('C3', 3), rep('C4', 5)),
                              term = c('negative regulation of mRNA processing',
                                       'alpha-beta T cell differentiation',
                                       'interleukin-4 production',
                                       'alpha-beta T cell activation',
                                       'cellular response to interleukin-4',
                                       'cytoplasmic translation',
                                       'RNA splicing',
                                       'ribosome biogenesis',
                                       'RNA localization',
                                       'cytoplasmic translation',
                                       'RNA transport',
                                       'RNA splicing'))

pdf('./age_negative_genes_time_cluster_heatmap.pdf', height = 8, width = 15)
visCluster(object = cm_filter, plot.type = "both", ht.col = c("#00c853", "white", "#c51162"),
           add.box = T, add.line = T, line.side = "left", show_row_dend = F,
           annoTerm.data = GOterm_selected, boxcol = ggsci::pal_npg()(3))
dev.off()

write.table(clini_nochi, './aging_cohortClini_nochi.txt', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(tpm_nochi, './aging_cohortTPM_nochi_forTIMER2.txt', col.names = T, row.names = T, quote = F, sep = '\t')

sampleTable <- clini_nochi
colnames(sampleTable)[2] <- 'condition'
sampleTable$condition <- gsub(' ', '_', sampleTable$condition)
sampleTable$condition <- factor(sampleTable$condition, levels = c('Youth', 'Middle_age', 'The_elderly'))

aging_counts <- read.table('./aging_cohortCounts.txt', header = T, sep = '\t', stringsAsFactors = F)
counts_nochi <- aging_counts[, clini_nochi$SampleID]
clini_nochi$SampleID == colnames(counts_nochi)
counts_nochi <- round(counts_nochi)

selected_sample <- c(1:16, 17:38)
dds <- DESeqDataSetFromMatrix(counts_nochi[, selected_sample], sampleTable[selected_sample, ], ~condition)
dds <- dds[rowSums(counts(dds)) >= ncol(counts_nochi[, selected_sample]), ]
vsd <- vst(dds)
dds <- DESeq(dds)
res <- results(dds)

Diff_Matrix <- res@listData
Diff_Matrix <- as.data.frame(Diff_Matrix)
row.names(Diff_Matrix) <- res@rownames

diff_noNA <- Diff_Matrix[rowSums(is.na(Diff_Matrix)) == 0, ]
diff_005 <- diff_noNA[diff_noNA$padj < 0.05, ]
diff_005_058 <- diff_005[abs(diff_005$log2FoldChange) > 0.58, ]
write.csv(diff_005_058, './MiddleAge_VS_Youth_diffgene_0.05_0.58.csv', col.names = T, row.names = T, quote = F)

high_exp <- rownames(diff_005_058)[diff_005_058$log2FoldChange > 0]
low_exp <- rownames(diff_005_058)[diff_005_058$log2FoldChange < 0]

GO_en_BP <- enrichGO(gene = high_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './MiddleAge_VS_Youth_highexp_GO_BP.csv', row.names = F, col.names = T, quote = F)

GO_en_BP <- enrichGO(gene = low_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './MiddleAge_VS_Youth_lowexp_GO_BP.csv', row.names = F, col.names = T, quote = F)

selected_sample <- c(17:38, 39:67)
dds <- DESeqDataSetFromMatrix(counts_nochi[, selected_sample], sampleTable[selected_sample, ], ~condition)
dds <- dds[rowSums(counts(dds)) >= ncol(counts_nochi[, selected_sample]), ]
vsd <- vst(dds)
dds <- DESeq(dds)
res <- results(dds)

Diff_Matrix <- res@listData
Diff_Matrix <- as.data.frame(Diff_Matrix)
row.names(Diff_Matrix) <- res@rownames

diff_noNA <- Diff_Matrix[rowSums(is.na(Diff_Matrix)) == 0, ]
diff_005 <- diff_noNA[diff_noNA$padj < 0.05, ]
diff_005_058 <- diff_005[abs(diff_005$log2FoldChange) > 0.58, ]
write.csv(diff_005_058, './Elderly_VS_MiddleAge_diffgene_0.05_0.58.csv', col.names = T, row.names = T, quote = F)

high_exp <- rownames(diff_005_058)[diff_005_058$log2FoldChange > 0]
low_exp <- rownames(diff_005_058)[diff_005_058$log2FoldChange < 0]

GO_en_BP <- enrichGO(gene = high_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './Elderly_VS_MiddleAge_highexp_GO_BP.csv', row.names = F, col.names = T, quote = F)

GO_en_BP <- enrichGO(gene = low_exp, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_en_BP_res <- GO_en_BP@result
write.csv(GO_en_BP_res, './Elderly_VS_MiddleAge_lowexp_GO_BP.csv', row.names = F, col.names = T, quote = F)

write.csv(diff_005_058, './MiddleAge_VS_Youth_diffgene_0.05_0.58.csv', col.names = T, row.names = T, quote = F)
write.csv(diff_005_058, './Elderly_VS_MiddleAge_diffgene_0.05_0.58.csv', col.names = T, row.names = T, quote = F)

middle_vs_youth <- read.csv('./MiddleAge_VS_Youth_diffgene_0.05_0.58.csv', header = T, stringsAsFactors = F)
elder_vs_middle <- read.csv('./Elderly_VS_MiddleAge_diffgene_0.05_0.58.csv', header = T, stringsAsFactors = F)

middle_vs_youth_high <- middle_vs_youth[middle_vs_youth$log2FoldChange > 0, ]
middle_vs_youth_low <- middle_vs_youth[middle_vs_youth$log2FoldChange < 0, ]

elder_vs_middle_high <- elder_vs_middle[elder_vs_middle$log2FoldChange > 0, ]
elder_vs_middle_low <- elder_vs_middle[elder_vs_middle$log2FoldChange < 0, ]

diff_bar_df <- data.frame(number = c(nrow(middle_vs_youth_high), nrow(middle_vs_youth_low), nrow(elder_vs_middle_high), nrow(elder_vs_middle_low)),
                          high_low = c('high', 'low', 'high', 'low'),
                          Group = c('MiddleAge_VS_Youth', 'MiddleAge_VS_Youth', 'Elderly_VS_MiddleAge', 'Elderly_VS_MiddleAge'))

pdf('./diffgene_3group_barplot_df.pdf', height = 5, width = 4)
ggplot(data = diff_bar_df, mapping = aes(x = Group, y = number, fill = high_low)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

senmayo <- read.csv('./2022_NC_cell_senescence_gene_set_GSVAscore_higher_cell_more_senescent_Human_PMID35974106.csv', header = T, stringsAsFactors = F)
senmayo_list <- list('SenMayo' = senmayo$Gene.human.)
saveRDS(senmayo_list, './SenMayo_list_forGSVA.rds')

cellage <- read.csv('./cellAge1.csv', header = T, sep = ';', stringsAsFactors = F)
cellage_induce <- cellage[cellage$senescence_effect == 'Induces', ]
cellage_inhibit <- cellage[cellage$senescence_effect == 'Inhibits', ]

cellage_induce_list <- list('CellAge_induceAging' = cellage_induce$gene_name)
cellage_inhibit_list <- list('CellAge_inhibitAging' = cellage_inhibit$gene_name)

saveRDS(cellage_induce_list, './CellAge_induceAging_list_forGSVA.rds')
saveRDS(cellage_inhibit_list, './CellAge_inhibitAging_list_forGSVA.rds')

immport_genes <- read.csv('./immuneGeneList_fromImmPort.csv', header = T, stringsAsFactors = F)
immport_final_list <- list()

immune_pathway <- immport_genes$Category
immune_pathway <- unique(immune_pathway)

for (i in 1:length(immune_pathway)) {
  pathway_name <- immune_pathway[i]
  pathway_gene_list <- immport_genes[immport_genes$Category == pathway_name, ]$Symbol
  immport_final_list[[pathway_name]] <- as.character(pathway_gene_list)
}

saveRDS(immport_final_list, './immPort_17_immunePathway_list_forGSVA.rds')

aging_immune_list <- c(senmayo_list, cellage_induce_list, cellage_inhibit_list, immport_final_list)

write.table(tpm_nochi, './aging_cohortTPM_nochi_forTIMER2.txt', col.names = T, row.names = T, quote = F, sep = '\t')

clini_nochi_ordered <- clini_nochi[order(clini_nochi$age), ]
tpm_nochi_ordered <- tpm_nochi[, clini_nochi_ordered$SampleID]
tpm_nochi_ordered_log2 <- log2(tpm_nochi_ordered + 1)
tpm_nochi_ordered_log2 <- as.matrix(tpm_nochi_ordered_log2)

aging_immune_gsva <- gsva(tpm_nochi_ordered_log2, aging_immune_list, verbose = TRUE, method = "ssgsea")

aging_immune_gsva <- as.data.frame(t(aging_immune_gsva))
age_col <- clini_nochi_ordered[3]

rownames(age_col) == rownames(aging_immune_gsva)

spearman_test <- corr.test(x = aging_immune_gsva, y = age_col, use = "pairwise", method = "spearman", adjust = "BH", ci = F)

cor_df <- as.data.frame(spearman_test$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(2, 1)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)

p_adj <- as.data.frame(spearman_test$p)
p_adj$ASE_id <- row.names(p_adj)
p_adj <- p_adj[, c(2, 1)]
p_adj_melt <- melt(p_adj, value.name = 'V2')
p_adj_melt <- unique(p_adj_melt)

colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
cor_melt$p_adjust <- p_adj_melt$V2
exp_age_spearman <- cor_melt
colnames(exp_age_spearman)[1:2] <- c('gene', 'age')
write.table(exp_age_spearman, './GSVAscore_VS_age_SpearmanResult.txt', col.names = T, row.names = F, sep = '\t', quote = F)

aging_immune_gsva_scatter <- aging_immune_gsva
aging_immune_gsva_scatter$age <- age_col$age

pdf('./Interferon_Receptor_VS_age_GSVA_ScatterPlot.pdf', height = 4, width = 5)
ggplot(aging_immune_gsva_scatter, aes(x = age, y = Interferon_Receptor)) + geom_point() +
  geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()

pdf('./CellAge_inhibitAging_VS_age_GSVA_ScatterPlot.pdf', height = 4, width = 5)
ggplot(aging_immune_gsva_scatter, aes(x = age, y = CellAge_inhibitAging)) + geom_point() +
  geom_smooth(method = "lm", color = "black", fill = "lightgray")
dev.off()