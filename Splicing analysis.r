library(psych)
library(reshape2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(cowplot)
library(gtable)
library(ggpubr)
library(grid)
library(UpSetR)
library(IsoformSwitchAnalyzeR)
library(tximport)
library(AnnotationDbi)

write.table(all_df_selected, './AgingCohort_noChi_PSImatrix.txt', col.names = T, row.names = T, sep = '\t', quote = F)

psi_mat <- read.table('./AgingCohort_noChi_PSImatrix.txt', header = T, sep = '\t', stringsAsFactors = F)
ase_num <- colSums(!is.na(psi_mat))
ase_num_df <- as.data.frame(ase_num)
ase_num_df$sample_id <- rownames(ase_num_df)
clin_info <- read.csv('./RNAgroup_aging_fromZhangQi.csv', header = T, stringsAsFactors = F)
clin_info <- clin_info[clin_info$SampleID %in% colnames(psi_mat), ]
rownames(clin_info) <- clin_info$SampleID
rownames(ase_num_df) == rownames(clin_info)
ase_num_df <- ase_num_df[rownames(clin_info), ]
rownames(ase_num_df) == rownames(clin_info)
ase_num_clin <- cbind(ase_num_df, clin_info)

median(ase_num_clin$ase_num[ase_num_clin$group == 'Youth'])
median(ase_num_clin$ase_num[ase_num_clin$group == 'Middle age'])
median(ase_num_clin$ase_num[ase_num_clin$group == 'The elderly'])
mean(ase_num_clin$ase_num[ase_num_clin$Gender == 'M'])
mean(ase_num_clin$ase_num[ase_num_clin$Gender == 'F'])
ase_num_clin <- ase_num_clin[, -2]
ase_num_clin$SampleID <- factor(ase_num_clin$SampleID, levels = clin_info$SampleID)

nrow(psi_mat)
psi_mat_type <- psi_mat
psi_mat_type$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', rownames(psi_mat_type))

psi_aa <- psi_mat_type[psi_mat_type$ASEtype == 'AA', ]
psi_ad <- psi_mat_type[psi_mat_type$ASEtype == 'AD', ]
psi_ap <- psi_mat_type[psi_mat_type$ASEtype == 'AP', ]
psi_at <- psi_mat_type[psi_mat_type$ASEtype == 'AT', ]
psi_es <- psi_mat_type[psi_mat_type$ASEtype == 'ES', ]
psi_me <- psi_mat_type[psi_mat_type$ASEtype == 'ME', ]
psi_ri <- psi_mat_type[psi_mat_type$ASEtype == 'RI', ]

num_aa <- colSums(!is.na(psi_aa[, -68]))
num_ad <- colSums(!is.na(psi_ad[, -68]))
num_ap <- colSums(!is.na(psi_ap[, -68]))
num_at <- colSums(!is.na(psi_at[, -68]))
num_es <- colSums(!is.na(psi_es[, -68]))
num_me <- colSums(!is.na(psi_me[, -68]))
num_ri <- colSums(!is.na(psi_ri[, -68]))

num_aa_df <- as.data.frame(num_aa)
num_ad_df <- as.data.frame(num_ad)
num_ap_df <- as.data.frame(num_ap)
num_at_df <- as.data.frame(num_at)
num_es_df <- as.data.frame(num_es)
num_me_df <- as.data.frame(num_me)
num_ri_df <- as.data.frame(num_ri)

num_aa_df$ASEtype <- 'AA'
num_ad_df$ASEtype <- 'AD'
num_ap_df$ASEtype <- 'AP'
num_at_df$ASEtype <- 'AT'
num_es_df$ASEtype <- 'ES'
num_me_df$ASEtype <- 'ME'
num_ri_df$ASEtype <- 'RI'

num_aa_df$SampleID <- rownames(num_aa_df)
num_ad_df$SampleID <- rownames(num_ad_df)
num_ap_df$SampleID <- rownames(num_ap_df)
num_at_df$SampleID <- rownames(num_at_df)
num_es_df$SampleID <- rownames(num_es_df)
num_me_df$SampleID <- rownames(num_me_df)
num_ri_df$SampleID <- rownames(num_ri_df)

colnames(num_aa_df)[1] <- 'Number'
colnames(num_ad_df)[1] <- 'Number'
colnames(num_ap_df)[1] <- 'Number'
colnames(num_at_df)[1] <- 'Number'
colnames(num_es_df)[1] <- 'Number'
colnames(num_me_df)[1] <- 'Number'
colnames(num_ri_df)[1] <- 'Number'

num_7type <- rbind(num_aa_df, num_ad_df, num_ap_df, num_at_df, num_es_df, num_me_df, num_ri_df)
num_7type$SampleID <- factor(num_7type$SampleID, levels = clin_info$SampleID)

seven_colors <- c('#c2a5cf', '#df65b0', '#ef6548', '#ffed6f', '#225ea8', '#cccc00', '#b8e186')

theme_set(theme_cowplot())
pdf('./AgingCohort_ASEProportion_7type.pdf', height = 4, width = 10)
ggplot(data = num_7type, mapping = aes(x = SampleID, y = Number, fill = ASEtype)) +
  geom_bar(stat = "identity", position = 'fill', width = 0.5, alpha = 1) +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  scale_fill_manual(values = seven_colors)
dev.off()

pdf('./AgingCohort_ASENumber_7type.pdf', height = 4, width = 10)
ggplot(data = num_7type, mapping = aes(x = SampleID, y = Number, fill = ASEtype)) +
  geom_bar(stat = "identity", position = 'stack', width = 0.5, alpha = 1) +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  scale_fill_manual(values = seven_colors)
dev.off()

table(psi_mat_type$ASEtype)
num_es_group <- num_es_df
num_es_group <- merge(num_es_group, clin_info, by = 'SampleID', all.x = T)
mean(num_es_group$Number[num_es_group$group == 'Youth'])
mean(num_es_group$Number[num_es_group$group == 'Middle age'])
mean(num_es_group$Number[num_es_group$group == 'The elderly'])

ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  hinvert_title_grob <- function(grob) {
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  index <- which(g2$layout$name == 'ylab-l')
  ylab <- g2$grobs[[index]]
  ylab <- hinvert_title_grob(ylab)
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  index <- which(g2$layout$name == 'axis-l')
  yaxis <- g2$grobs[[index]]
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  yaxis$children[[2]] <- ticks
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

ase_num_clin <- ase_num_clin[order(ase_num_clin$age), ]
ase_num_clin$SampleID <- factor(ase_num_clin$SampleID, levels = ase_num_clin$SampleID)
num_7type$SampleID <- factor(num_7type$SampleID, levels = ase_num_clin$SampleID)
num_7type$ASEtype <- factor(num_7type$ASEtype, levels = rev(c('AA', 'AD', 'AP', 'AT', 'ES', 'ME', 'RI')))
levels(ase_num_clin$SampleID) == levels(num_7type$SampleID)

theme_set(theme_cowplot())
pdf("./ASE_numberInEachAgingCohortSample.pdf", height = 3, width = 5)
p1 <- ggplot(data = num_7type, mapping = aes(x = SampleID, y = Number, fill = ASEtype)) +
  geom_bar(stat = "identity", position = 'fill', width = 0.5, alpha = 1) +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  scale_fill_manual(values = rev(seven_colors))
p2 <- ggplot(data = ase_num_clin, mapping = aes(x = SampleID, y = ase_num)) + geom_point(size = 0.5)
ggplot2.two_y_axis(p1, p2)
dev.off()

cor.test(ase_num_clin$ase_num, ase_num_clin$age, method = 'spearman')

write.table(all_df_selected, './AgingCohort_noChi_PSImatrix.txt', col.names = T, row.names = T, sep = '\t', quote = F)

psi_mat <- read.table('./AgingCohort_noChi_PSImatrix.txt', header = T, sep = '\t', stringsAsFactors = F)
nrow(psi_mat)
clin_info <- read.csv('./RNAgroup_aging_fromZhangQi.csv', header = T, stringsAsFactors = F)
clin_info <- clin_info[clin_info$SampleID %in% colnames(psi_mat), ]
rownames(clin_info) <- clin_info$SampleID
psi_mat <- psi_mat[, clin_info$SampleID]
colnames(psi_mat) == clin_info$SampleID

psi_filt <- psi_mat[(rowSums(!is.na(psi_mat[, 1:16])) >= 2 &
                     rowSums(!is.na(psi_mat[, 17:38])) >= 2 &
                     rowSums(!is.na(psi_mat[, 39:67])) >= 2), ]

exp_mat <- t(psi_filt)
age_df <- clin_info[, c(1, 3)]
age_df <- age_df[rownames(exp_mat), ]
all(rownames(exp_mat) == rownames(age_df))
age_df <- age_df[2]
exp_mat <- as.matrix(exp_mat)
age_df <- as.matrix(age_df)

spearman_test <- corr.test(x = exp_mat, y = age_df, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df <- as.data.frame(spearman_test$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(2, 1)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)
p_df <- as.data.frame(spearman_test$p)
p_df$ASE_id <- row.names(p_df)
p_df <- p_df[, c(2, 1)]
p_melt <- melt(p_df, value.name = 'V2')
p_melt <- unique(p_melt)
colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
cor_melt$p_adjust <- p_melt$V2
exp_age_spearman <- cor_melt
colnames(exp_age_spearman)[1:2] <- c('gene', 'age')
colnames(psi_filt) == clin_info$SampleID
diff_mean <- rowMeans(psi_filt[, 39:67], na.rm = T) - rowMeans(psi_filt[, 1:16], na.rm = T)
diff_mean_df <- as.data.frame(diff_mean)
all(exp_age_spearman$gene == rownames(diff_mean_df))
exp_age_spearman$DiffMean_OldMinusYoung <- diff_mean_df$diff_mean
exp_age_spearman_filtered <- exp_age_spearman[(exp_age_spearman$p_adjust < 0.05 & abs(exp_age_spearman$DiffMean_OldMinusYoung) > 0.1), ]

age_ase <- exp_age_spearman_filtered
age_ase$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', age_ase$gene)

aa_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'AA'])
ad_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'AD'])
ap_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'AP'])
at_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'AT'])
es_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'ES'])
me_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'ME'])
ri_genes <- unique(age_ase$relatedGene[age_ase$ASEtype == 'RI'])

aa_num <- length(aa_genes)
ad_num <- length(ad_genes)
ap_num <- length(ap_genes)
at_num <- length(at_genes)
es_num <- length(es_genes)
me_num <- length(me_genes)
ri_num <- length(ri_genes)

gene_df <- data.frame(Var1 = c('AA', 'AD', 'AP', 'AT', 'ES', 'ME', 'RI'),
                      Freq = c(aa_num, ad_num, ap_num, at_num, es_num, me_num, ri_num))
gene_df$Variables <- 'Gene'
ase_df <- as.data.frame(table(age_ase$ASEtype))
ase_df$Variables <- 'ASE'
ase_gene_df <- rbind(ase_df, gene_df)
colnames(ase_gene_df) <- c('type', 'Value', 'Variables')

theme_set(theme_cowplot())
pdf("./aging_associated_ASE_ASEtype_GeneNumber.pdf", height = 3, width = 5)
ggplot(data = ase_gene_df, aes(x = type, y = Value, fill = Variables)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8) +
  labs(x = "ASE Type", y = "Numbers") + theme(legend.key = element_blank()) +
  geom_text(aes(label = paste(Value)), angle = 15, vjust = -0.3, position = position_dodge(0.9), size = 2.8) +
  scale_fill_manual(values = c("steelblue3", "tomato3"))
dev.off()

list_for_upset <- list(AA = aa_genes, AD = ad_genes, AP = ap_genes, AT = at_genes, ES = es_genes, ME = me_genes, RI = ri_genes)

pdf(file = "./UpsetPlot.pdf", width = 10, height = 6)
upset(fromList(list_for_upset), nsets = 7, number.angles = 15, point.size = 3, line.size = 0.7,
      mainbar.y.label = "Gene Intersections", sets.x.label = "Set Size", mb.ratio = c(0.6, 0.4),
      text.scale = c(1.8, 1.3, 1.8, 1.3, 1.5, 1.5), empty.intersections = "off")
dev.off()

age_mrna <- read.table('./age_associated_ProteinCodingGene_noChi_Gencode_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)

data(splicingFactors)
sum(age_mrna$gene %in% splicingFactors$GeneSymbol)
splicingFactors$GeneSymbol <- gsub('SFRS', 'SRSF', splicingFactors$GeneSymbol)
age_sf <- age_mrna[age_mrna$gene %in% splicingFactors$GeneSymbol, ]

aging_cohortTPM <- read.table('./aging_cohortTPM.txt', header = T, sep = '\t', stringsAsFactors = F)
aging_cohortTPM_sf <- aging_cohortTPM[age_sf$gene, ]
aging_cohortTPM_sf_nochi <- aging_cohortTPM_sf[, grep('Chi', invert = T, colnames(aging_cohortTPM_sf))]
psi_mat_nochi <- read.table('./AgingCohort_noChi_PSImatrix.txt', header = T, sep = '\t', stringsAsFactors = F)
age_ase_psi <- psi_mat_nochi[age_ase$gene, ]
age_ase_psi <- age_ase_psi[, colnames(aging_cohortTPM_sf_nochi)]
colnames(aging_cohortTPM_sf_nochi) == colnames(age_ase_psi)
aging_cohortTPM_sf_nochi_t <- t(aging_cohortTPM_sf_nochi)
age_ase_psi_t <- t(age_ase_psi)

spearman_test <- corr.test(x = age_ase_psi_t, y = aging_cohortTPM_sf_nochi_t, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df <- as.data.frame(spearman_test$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(53, 1:52)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)
p_df <- as.data.frame(spearman_test$p)
p_df$ASE_id <- row.names(p_df)
p_df <- p_df[, c(53, 1:52)]
p_melt <- melt(p_df, value.name = 'V2')
p_melt <- unique(p_melt)
colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
all(cor_melt$diffASE_id == p_melt$ASE_id)
cor_melt$p_adjust <- p_melt$V2
rbp_ase_spearman <- cor_melt
rbp_ase_spearman_filt <- rbp_ase_spearman[(abs(rbp_ase_spearman$correlation) > 0.5 & rbp_ase_spearman$p_adjust < 0.05), ]
length(unique(rbp_ase_spearman_filt$diffRBP))
length(unique(rbp_ase_spearman_filt$diffASE_id))

diffASE_id <- unique(rbp_ase_spearman_filt[1])
diffASE_type_gene <- data.frame(diffASE = diffASE_id)
diffASE_type_gene$Gene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene$diffASE_id)
diffASE_type_gene$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene$diffASE_id)
colnames(diffASE_type_gene) <- c('diffASE', 'Gene', 'ASEtype')
diffASE_related_gene <- diffASE_type_gene[2]
diffASE_type_gene <- diffASE_type_gene[, c(1, 3)]
rbp_ase_spearman_filt <- merge(rbp_ase_spearman_filt, diffASE_type_gene, by.x = 'diffASE_id', by.y = 'diffASE', all.x = T)
rbp_ase_spearman_filt$regulate_direct <- ifelse(rbp_ase_spearman_filt$correlation > 0, 'Positive', 'Negative')
diffASE_related_gene <- unique(diffASE_related_gene)
length(diffASE_related_gene$Gene)

go_en <- enrichGO(gene = diffASE_related_gene$Gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_en_result <- go_en@result

down_ase <- age_sf[age_sf$Direction == 'downregulated', ]
down_ase_genes <- down_ase$gene
down_rbp_ase_filt <- rbp_ase_spearman_filt[rbp_ase_spearman_filt$diffRBP %in% down_ase_genes, ]
diffASE_id_down <- unique(down_rbp_ase_filt[1])
diffASE_type_gene_down <- data.frame(diffASE = diffASE_id_down)
diffASE_type_gene_down$Gene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene_down$diffASE_id)
diffASE_type_gene_down$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene_down$diffASE_id)
colnames(diffASE_type_gene_down) <- c('diffASE', 'Gene', 'ASEtype')
diffASE_related_gene_down <- diffASE_type_gene_down[2]
diffASE_type_gene_down <- diffASE_type_gene_down[, c(1, 3)]
diffASE_related_gene_down <- unique(diffASE_related_gene_down)
length(diffASE_related_gene_down$Gene)

go_en_down <- enrichGO(gene = diffASE_related_gene_down$Gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_en_result_down <- go_en_down@result

cyto_trans_genes <- go_en_result_down$geneID[go_en_result_down$Description == 'cytoplasmic translation']
cyto_trans_genes <- strsplit(cyto_trans_genes, split = '/')[[1]]

rbp_ase_spearman_filt$diffASE_relatedGene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', rbp_ase_spearman_filt$diffASE_id)
sf_vs_cyto <- rbp_ase_spearman_filt[rbp_ase_spearman_filt$diffASE_relatedGene %in% cyto_trans_genes, ]
length(unique(sf_vs_cyto$diffRBP))
length(unique(sf_vs_cyto$diffASE_relatedGene))

age_ase$relatedGene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', age_ase$gene)
age_ase_related_genes <- unique(age_ase$relatedGene)

go_en_ase <- enrichGO(gene = age_ase_related_genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_en_result_ase <- go_en_ase@result

selected_desc <- c('protein modification by small protein removal', 'ribonucleoprotein complex assembly', 'protein deubiquitination', 'cytoplasmic translation')
go_en_result_selected <- go_en_result_ase[go_en_result_ase$Description %in% selected_desc, ]
go_en_result_selected$minus_log10_pvalue <- -log10(go_en_result_selected$pvalue)
go_en_result_selected$Description <- factor(go_en_result_selected$Description, levels = rev(go_en_result_selected$Description))

theme_set(theme_cowplot())
pdf('./age_associated_ASE_relatedGene_GOresult_barplot.pdf', height = 4, width = 8)
ggplot(go_en_result_selected, aes(x = minus_log10_pvalue, y = Description)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8)
dev.off()

aging_base_dir <- './salmon_output'
aging_files <- list.files(aging_base_dir)
aging_files <- aging_files[grep('[Z|Y]', aging_files)]
aging_file_list <- paste0(aging_base_dir, '/', aging_files, '/quant.sf')
all_file_list <- aging_file_list
all_file_list_names <- gsub('.*\\/(.*)/quant.sf', '\\1', all_file_list)
names(all_file_list) <- all_file_list_names
tx2gene <- read.table('./tx2gene.gencode_GRch38_v30.txt', header = T, sep = '\t', stringsAsFactors = F)

all_transcript_exp <- tximport(files = all_file_list, type = "salmon", txOut = T, txIn = T, countsFromAbundance = 'dtuScaledTPM', tx2gene = tx2gene)
tx_tpm <- all_transcript_exp$abundance
rownames(tx_tpm) <- gsub('(.*)\\..*', '\\1', rownames(tx_tpm))
tx_counts <- all_transcript_exp$counts
rownames(tx_counts) <- gsub('(.*)\\..*', '\\1', rownames(tx_counts))

tx2gene$gene_id_clean <- gsub('(.*)\\..*', '\\1', tx2gene$gene_id)
tx2gene$txid_clean <- gsub('(.*)\\..*', '\\1', tx2gene$txid)
tx2gene$gene_symbol <- mapIds(org.Hs.eg.db, keys = tx2gene$gene_id_clean, keytype = 'ENSEMBL', column = 'SYMBOL')
age_ase_related_gene_tx <- tx2gene[tx2gene$gene_symbol %in% age_ase_related_genes, ]
selected_txid <- age_ase_related_gene_tx$txid_clean[age_ase_related_gene_tx$txid_clean %in% rownames(tx_tpm)]
related_gene_tx_tpm <- tx_tpm[selected_txid, ]
related_gene_tx_tpm <- as.data.frame(related_gene_tx_tpm)
related_gene_tx_counts <- tx_counts[selected_txid, ]
related_gene_tx_counts <- as.data.frame(related_gene_tx_counts)
related_gene_tx_tpm$isoform_id <- rownames(related_gene_tx_tpm)
related_gene_tx_counts$isoform_id <- rownames(related_gene_tx_counts)
related_gene_tx_tpm <- related_gene_tx_tpm[, c(68, 1:67)]
related_gene_tx_counts <- related_gene_tx_counts[, c(68, 1:67)]
rownames(related_gene_tx_tpm) <- 1:nrow(related_gene_tx_tpm)
rownames(related_gene_tx_counts) <- 1:nrow(related_gene_tx_counts)

clin_info_nochi <- read.table('./aging_cohortCliniInfo_noChi.txt', header = T, sep = '\t', stringsAsFactors = F)
clin_info_nochi <- clin_info_nochi[order(clin_info_nochi$age), ]
clin_info_nochi <- clin_info_nochi[, c(1, 2)]
clin_info_nochi$group[clin_info_nochi$group == 'Middle age'] <- 'The elderly'
design_mat <- clin_info_nochi
colnames(design_mat) <- c('sampleID', 'condition')

gtf <- read.table("./gencode_GRCh38_v30.gtf", skip = 5, header = F, sep = '\t', stringsAsFactors = F)
gtf_sep <- separate(gtf, col = V9, into = paste0('al', 1:2), sep = ';')
gtf_sep$V9 <- gtf$V9
gtf_sep_again <- separate(gtf_sep, col = al2, into = paste0('b1', 1:5), sep = ' ')
tx_id_with_dot <- gtf_sep_again$b13
tx_id_with_dot <- unique(tx_id_with_dot)
tx_id_with_dot <- tx_id_with_dot[grep('ENST', tx_id_with_dot)]
tx_id_no_dot <- gsub('(.*)\\..*', '\\1', tx_id_with_dot)
tx_id_dot_nodot_df <- data.frame(dot = tx_id_with_dot, nodot = tx_id_no_dot)
selected_txid_df <- data.frame(nodot = selected_txid)
selected_txid_df <- merge(selected_txid_df, tx_id_dot_nodot_df, by = 'nodot', all.x = T)
rownames(selected_txid_df) <- selected_txid_df$nodot
selected_txid_df <- selected_txid_df[related_gene_tx_tpm$isoform_id, ]
all(selected_txid_df$nodot == related_gene_tx_tpm$isoform_id)
related_gene_tx_tpm$isoform_id <- selected_txid_df$dot
all(selected_txid_df$nodot == related_gene_tx_counts$isoform_id)
related_gene_tx_counts$isoform_id <- selected_txid_df$dot

gtf_sep_again_filt <- gtf_sep_again[gtf_sep_again$b13 %in% selected_txid_df$dot, ]
final_gtf <- gtf_sep_again_filt[, paste0('V', 1:9)]

fasta <- "./gencode.v30.transcripts.fa"

aSwitchList <- importRdata(
  isoformCountMatrix   = related_gene_tx_counts,
  isoformRepExpression = related_gene_tx_tpm,
  designMatrix         = design_mat,
  isoformExonAnnoation = final_gtf,
  isoformNtFasta = fasta,
  showProgress = FALSE
)

saveRDS(aSwitchList, './ASErelatedGene_aSwitchList.rds')
switchList_part1 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  pathToOutput = "./",
  outputSequences      = TRUE,
  prepareForWebServers = FALSE
)

down_sf <- c('PTBP1', 'HNRNPM', 'HNRNPK', 'HNRNPA2B1', 'U2AF2', 'PRPF19', 'RNPS1', 'SF1', 'SNRPA1', 'SF3B2', 'PCBP1')

aging_cohortTPM <- read.table('./aging_cohortTPM.txt', header = T, sep = '\t', stringsAsFactors = F)
aging_cohortTPM_sf <- aging_cohortTPM[down_sf, ]
aging_cohortTPM_sf_nochi <- aging_cohortTPM_sf[, grep('Chi', invert = T, colnames(aging_cohortTPM_sf))]
psi_mat_nochi <- read.table('./AgingCohort_noChi_PSImatrix.txt', header = T, sep = '\t', stringsAsFactors = F)
age_ase_psi <- psi_mat_nochi[age_ase$gene, ]
age_ase_psi <- age_ase_psi[, colnames(aging_cohortTPM_sf_nochi)]
colnames(aging_cohortTPM_sf_nochi) == colnames(age_ase_psi)
aging_cohortTPM_sf_nochi_t <- t(aging_cohortTPM_sf_nochi)
age_ase_psi_t <- t(age_ase_psi)

spearman_test <- corr.test(x = age_ase_psi_t, y = aging_cohortTPM_sf_nochi_t, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df <- as.data.frame(spearman_test$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(12, 1:11)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)
p_df <- as.data.frame(spearman_test$p)
p_df$ASE_id <- row.names(p_df)
p_df <- p_df[, c(12, 1:11)]
p_melt <- melt(p_df, value.name = 'V2')
p_melt <- unique(p_melt)
colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
all(cor_melt$diffASE_id == p_melt$ASE_id)
cor_melt$p_adjust <- p_melt$V2
rbp_ase_spearman <- cor_melt
rbp_ase_spearman_filt <- rbp_ase_spearman[(abs(rbp_ase_spearman$correlation) > 0.5 & rbp_ase_spearman$p_adjust < 0.05), ]
table(rbp_ase_spearman_filt$diffRBP)

rbp_ase_spearman_filt_plot <- rbp_ase_spearman_filt
rbp_ase_spearman_filt_plot$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', rbp_ase_spearman_filt_plot$diffASE_id)
df_plot <- as.data.frame(table(rbp_ase_spearman_filt_plot$diffRBP, rbp_ase_spearman_filt_plot$ASEtype))
num_df <- as.data.frame(table(rbp_ase_spearman_filt_plot$diffRBP))
num_df <- num_df[order(num_df$Freq, decreasing = T), ]
levels_sf <- num_df$Var1
df_plot$Var1 <- factor(df_plot$Var1, levels = levels_sf)

theme_set(theme_cowplot())
pdf('./Down11SF_regulatingASE_number.pdf', height = 4, width = 5.5)
ggplot(data = df_plot, mapping = aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = 'stack', width = 0.7, alpha = 1) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  scale_fill_manual(values = seven_colors)
dev.off()

diffASE_id <- unique(rbp_ase_spearman_filt[1])
diffASE_type_gene <- data.frame(diffASE = diffASE_id)
diffASE_type_gene$Gene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene$diffASE_id)
diffASE_type_gene$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', diffASE_type_gene$diffASE_id)
colnames(diffASE_type_gene) <- c('diffASE', 'Gene', 'ASEtype')
diffASE_related_gene <- diffASE_type_gene[2]
diffASE_type_gene <- diffASE_type_gene[, c(1, 3)]
rbp_ase_spearman_filt <- merge(rbp_ase_spearman_filt, diffASE_type_gene, by.x = 'diffASE_id', by.y = 'diffASE', all.x = T)
rbp_ase_spearman_filt$regulate_direct <- ifelse(rbp_ase_spearman_filt$correlation > 0, 'Positive', 'Negative')
diffASE_related_gene <- unique(diffASE_related_gene)
length(diffASE_related_gene$Gene)

go_en <- enrichGO(gene = diffASE_related_gene$Gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_en_result <- go_en@result

age_ase$relatedGene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', age_ase$gene)
age_ase_related_genes <- unique(age_ase$relatedGene)

go_en_ase <- enrichGO(gene = age_ase_related_genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
go_en_result_ase <- go_en_ase@result

protein_mod_genes <- go_en_result_ase$geneID[go_en_result_ase$Description == 'protein modification by small protein removal']
ribo_assembly_genes <- go_en_result_ase$geneID[go_en_result_ase$Description == 'ribonucleoprotein complex assembly']
protein_deubiq_genes <- go_en_result_ase$geneID[go_en_result_ase$Description == 'protein deubiquitination']
cyto_trans_genes <- go_en_result_ase$geneID[go_en_result_ase$Description == 'cytoplasmic translation']
protein_mod_genes <- strsplit(protein_mod_genes, split = '/')[[1]]
ribo_assembly_genes <- strsplit(ribo_assembly_genes, split = '/')[[1]]
protein_deubiq_genes <- strsplit(protein_deubiq_genes, split = '/')[[1]]
cyto_trans_genes <- strsplit(cyto_trans_genes, split = '/')[[1]]
all_protein_genes <- c(protein_mod_genes, ribo_assembly_genes, protein_deubiq_genes, cyto_trans_genes)
all_protein_genes <- unique(all_protein_genes)

psi_protein <- psi_filt
psi_protein$gene <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', rownames(psi_protein))
psi_protein <- psi_protein[psi_protein$gene %in% all_protein_genes, ]
psi_protein <- psi_protein[, -ncol(psi_protein)]

clin_info_prot <- read.table('./Proteome_CliniInfo_nochi.txt', header = T, sep = '\t', stringsAsFactors = F)
prot_norm <- read.csv('./ZD_Lumos_ZD018_All_65Sample_Merge_Protein_UP1_Normalization.csv', header = T, stringsAsFactors = F)
prot_norm_nochi <- prot_norm[, grep('XC', invert = T, colnames(prot_norm))]
prot_norm_nochi <- prot_norm_nochi[rowSums(!is.na(prot_norm_nochi[, 4:51])) > 0, ]
prot_norm_nochi <- prot_norm_nochi[!duplicated(prot_norm_nochi$Symbol), ]
rownames(prot_norm_nochi) <- prot_norm_nochi$Symbol
prot_norm_nochi <- prot_norm_nochi[, -c(1, 2, 3)]
all(clin_info_prot$ID %in% colnames(prot_norm_nochi))
prot_norm_nochi <- prot_norm_nochi[, clin_info_prot$ID]

inter_samples <- colnames(prot_norm_nochi)[colnames(prot_norm_nochi) %in% colnames(psi_protein)]
prot_norm_selected <- prot_norm_nochi[, inter_samples]
psi_protein_selected <- psi_protein[, inter_samples]

prot_norm_selected <- as.matrix(t(prot_norm_selected))
psi_protein_selected <- as.matrix(t(psi_protein_selected))
rownames(psi_protein_selected) == rownames(prot_norm_selected)

spearman_test <- corr.test(x = psi_protein_selected, y = prot_norm_selected, use = "pairwise", method = "spearman", adjust = "BH", ci = F)
cor_df <- as.data.frame(spearman_test$r)
cor_df$ASE_id <- row.names(cor_df)
n_col <- ncol(cor_df)
cor_df <- cor_df[, c(n_col, 1:(n_col - 1))]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)
p_df <- as.data.frame(spearman_test$p)
p_df$ASE_id <- row.names(p_df)
p_df <- p_df[, c(n_col, 1:(n_col - 1))]
p_melt <- melt(p_df, value.name = 'V2')
p_melt <- unique(p_melt)
colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
cor_melt$p_adjust <- p_melt$V2
exp_age_spearman <- cor_melt
colnames(exp_age_spearman)[1:2] <- c('gene', 'age')
exp_age_spearman_filt <- exp_age_spearman[exp_age_spearman$p_adjust < 0.05, ]
exp_age_spearman_filt <- exp_age_spearman_filt[!is.na(exp_age_spearman_filt$correlation), ]
unique(exp_age_spearman_filt$gene)
length(unique(exp_age_spearman_filt$age))

exp_age_spearman_filt_07 <- exp_age_spearman_filt[abs(exp_age_spearman_filt$correlation) > 0.7, ]
unique(exp_age_spearman_filt$gene)
length(unique(exp_age_spearman_filt_07$age))
exp_age_spearman_filt_07$correlation_group <- 'positive'
exp_age_spearman_filt_07$correlation_group[exp_age_spearman_filt_07$correlation < 0] <- 'negative'
table(exp_age_spearman_filt_07$correlation_group)

exp_age_spearman_filt_07$ASEsymbol <- gsub('(.*)\\_.*\\_.*\\_.*\\_.*', '\\1', exp_age_spearman_filt_07$gene)
exp_age_spearman_filt_07$ASEtype <- gsub('.*\\_(.*)\\_.*\\_.*\\_.*', '\\1', exp_age_spearman_filt_07$gene)

protein_mod_df <- data.frame(ID = protein_mod_genes, GOterm = 'protein_modification_gene')
ribo_assembly_df <- data.frame(ID = ribo_assembly_genes, GOterm = 'ribonucleoprotein_complex_assembly_gene')
protein_deubiq_df <- data.frame(ID = protein_deubiq_genes, GOterm = 'protein_deubiquitination_gene')
cyto_trans_df <- data.frame(ID = cyto_trans_genes, GOterm = 'cytoplasmic_translation_gene')
ase_protein_df <- rbind(protein_mod_df, ribo_assembly_df, protein_deubiq_df, cyto_trans_df)
ase_protein_df$GOterm <- gsub('_gene', '', ase_protein_df$GOterm)
selected_freq <- as.data.frame(table(ase_protein_df$ID))
only_one <- selected_freq$Var1[selected_freq$Freq == 1]
only_one <- unique(only_one)
only_two <- selected_freq$Var1[selected_freq$Freq == 2]
only_two <- unique(only_two)

ase_protein_df_one <- ase_protein_df[ase_protein_df$ID %in% only_one, ]
ase_protein_df_two <- ase_protein_df[ase_protein_df$ID %in% only_two, ]
ase_protein_df_two <- ase_protein_df_two[order(ase_protein_df_two$ID), ]
rownames(ase_protein_df_two) <- 1:nrow(ase_protein_df_two)
for (i in seq(from = 1, to = nrow(ase_protein_df_two), by = 2)) {
  ase_protein_df_two$GOterm[i] <- paste(ase_protein_df_two$GOterm[i], ase_protein_df_two$GOterm[i + 1], sep = '; ')
}
ase_protein_df_two <- ase_protein_df_two[seq(from = 1, to = nrow(ase_protein_df_two), by = 2), ]
ase_protein_df_one_two <- rbind(ase_protein_df_one, ase_protein_df_two)
rownames(ase_protein_df_one_two) <- 1:nrow(ase_protein_df_one_two)
exp_age_spearman_filt_07 <- merge(exp_age_spearman_filt_07, ase_protein_df_one_two, by.x = 'ASEsymbol', by.y = 'ID', all.x = T)
df_freq <- as.data.frame(table(exp_age_spearman_filt_07$gene))
exp_age_spearman_filt_07 <- merge(exp_age_spearman_filt_07, df_freq, by.x = 'gene', by.y = 'Var1', all.x = T)

exp_age_spearman_filt_07_remove1 <- exp_age_spearman_filt_07
exp_age_spearman_filt_07_remove1$correlation <- round(exp_age_spearman_filt_07_remove1$correlation, digits = 2)
exp_age_spearman_filt_07_remove1$p_adjust <- round(exp_age_spearman_filt_07_remove1$p_adjust, digits = 2)
exp_age_spearman_filt_07_remove1 <- exp_age_spearman_filt_07_remove1[!(abs(exp_age_spearman_filt_07_remove1$correlation) == 1), ]
exp_age_spearman_filt_07_remove1 <- exp_age_spearman_filt_07_remove1[, -9]
df_freq_new <- as.data.frame(table(exp_age_spearman_filt_07_remove1$gene))
length(unique(exp_age_spearman_filt_07_remove1$age))
table(exp_age_spearman_filt_07_remove1$ASEtype)
7301 / sum(as.data.frame(table(exp_age_spearman_filt_07_remove1$ASEtype))$Freq)
exp_age_spearman_filt_07_remove1 <- merge(exp_age_spearman_filt_07_remove1, df_freq_new, by.x = 'gene', by.y = 'Var1', all.x = T)