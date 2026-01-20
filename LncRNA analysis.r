library(reshape2)
library(tidyr)
library(psych)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)

lncipedia <- read.delim('./lncipedia_version5.2_HighConfidenceSet_hg38.gtf', header = F, sep = '\t', stringsAsFactors = F, skip = 2)

lncipedia_id <- lncipedia[9]
lncipedia_id <- unique(lncipedia_id)
rownames(lncipedia_id) <- 1:nrow(lncipedia_id)
lncipedia_chr <- lncipedia[, 1:8]

into <- paste0("A", 2:150)
lncipedia_sep <- separate(lncipedia_id, col = V9, into = into, sep = ";")
lncipedia_sep <- as.matrix(lncipedia_sep)
dim(lncipedia_sep) <- c(nrow(lncipedia_sep) * ncol(lncipedia_sep), 1)
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- lncipedia_sep[!is.na(lncipedia_sep$V1), ]
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- unique(lncipedia_sep)
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- lncipedia_sep[grep('gene', lncipedia_sep$lncipedia_sep), ]
lncipedia_sep <- as.data.frame(lncipedia_sep)

into <- paste0("A", 1:50)
lncipedia_sep <- separate(lncipedia_sep, col = lncipedia_sep, into = into, sep = "\t")
lncipedia_sep <- as.matrix(lncipedia_sep)
dim(lncipedia_sep) <- c(nrow(lncipedia_sep) * ncol(lncipedia_sep), 1)
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- lncipedia_sep[!is.na(lncipedia_sep$V1), ]
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- lncipedia_sep[grep('gene', lncipedia_sep$lncipedia_sep), ]
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- unique(lncipedia_sep)
lncipedia_geneid <- lncipedia_sep[grep('gene_id', lncipedia_sep$lncipedia_sep), ]
lncipedia_alias <- lncipedia_sep[grep('gene_alias', lncipedia_sep$lncipedia_sep), ]
(length(lncipedia_geneid) + length(lncipedia_alias)) == nrow(lncipedia_sep)

lncipedia_geneid <- as.data.frame(lncipedia_geneid)
lncipedia_alias <- as.data.frame(lncipedia_alias)
lncipedia_geneid$lncipedia_geneid <- gsub('gene_id ', '', lncipedia_geneid$lncipedia_geneid)

into <- paste0("A", 1:5)
lncipedia_alias <- separate(lncipedia_alias, col = lncipedia_alias, into = into, sep = " ")
lncipedia_alias <- lncipedia_alias[3]
colnames(lncipedia_alias) <- 'lncipedia_geneid'

lnc_list <- rbind(lncipedia_geneid, lncipedia_alias)
lnc_list <- unique(lnc_list)

'LINC00324' %in% lnc_list$lncipedia_geneid

sum(lnc_list$lncipedia_geneid == 'LINC00324')

lnc_list$lncipedia_geneid <- gsub(' ', '', lnc_list$lncipedia_geneid)

coding_genes <- read.table('./ProteinCodingGene_list_from_gencodeV41.txt', header = T, sep = '\t', stringsAsFactors = F)

sum(lnc_list$lncipedia_geneid %in% coding_genes$gene_name)

lnc_list$lncipedia_geneid[lnc_list$lncipedia_geneid %in% coding_genes$gene_name]

lnc_list <- lnc_list[!(lnc_list$lncipedia_geneid %in% coding_genes$gene_name), ]

lnc_list <- as.data.frame(lnc_list)

lncipedia <- read.delim('./lncipedia_version5.2_HighConfidenceSet_hg38.gtf', header = F, sep = '\t', stringsAsFactors = F, skip = 2)
lncipedia_id <- lncipedia[9]
into <- paste0("A", 2:150)
lncipedia_sep <- separate(lncipedia_id, col = V9, into = into, sep = ";")
lncipedia_sep <- as.matrix(lncipedia_sep)
dim(lncipedia_sep) <- c(nrow(lncipedia_sep) * ncol(lncipedia_sep), 1)
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- lncipedia_sep[grep('gene_id', lncipedia_sep$V1), ]
lncipedia_sep <- as.data.frame(lncipedia_sep)
lncipedia_sep <- unique(lncipedia_sep)
lncipedia_sep$lncipedia_sep <- gsub('gene_id ', '', lncipedia_sep$lncipedia_sep)
lncipedia_sep$lncipedia_sep <- gsub(' ', '', lncipedia_sep$lncipedia_sep)
colnames(lncipedia_sep) <- 'LncRNA_geneID'

age_lnc <- read.table('./age_associated_lncRNA_noChi_LNCipedia_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)
age_coding <- read.table('./age_associated_ProteinCodingGene_noChi_Gencode_SpearmanResult.txt', header = T, sep = '\t', stringsAsFactors = F)

aging_tpm <- read.table('./aging_cohortTPM.txt', header = T, sep = '\t', stringsAsFactors = F)
tpm_nochi <- aging_tpm[, grep('Chi', invert = T, colnames(aging_tpm))]

tpm_lnc <- tpm_nochi[age_lnc$gene, ]
tpm_coding <- tpm_nochi[age_coding$gene, ]

all(colnames(tpm_coding) == colnames(tpm_lnc))

tpm_lnc <- t(tpm_lnc)
tpm_coding <- t(tpm_coding)

spearman_res <- corr.test(x = tpm_lnc, y = tpm_coding, use = "pairwise", method = "spearman", adjust = "BH", ci = F)

cor_df <- as.data.frame(spearman_res$r)
cor_df$ASE_id <- row.names(cor_df)
cor_df <- cor_df[, c(2210, 1:2209)]
cor_melt <- melt(cor_df, value.name = 'V2')
cor_melt <- unique(cor_melt)

p_adj <- as.data.frame(spearman_res$p)
p_adj$ASE_id <- row.names(p_adj)
p_adj <- p_adj[, c(2210, 1:2209)]
p_adj_melt <- melt(p_adj, value.name = 'V2')
p_adj_melt <- unique(p_adj_melt)

colnames(cor_melt) <- c('diffASE_id', 'diffRBP', 'correlation')
cor_melt$p_adjust <- p_adj_melt$V2
lnc_mRNA_spearman <- cor_melt
colnames(lnc_mRNA_spearman)[1:2] <- c('lncRNA', 'mRNA')

up_lnc <- age_lnc[age_lnc$correlation > 0, ]
up_lnc$Direction <- 'Upregulated'
down_lnc <- age_lnc[age_lnc$correlation < 0, ]
down_lnc$Direction <- 'Downregulated'

lnc_mRNA_spearman_090 <- lnc_mRNA_spearman[(abs(lnc_mRNA_spearman$correlation) >= 0.9 & lnc_mRNA_spearman$p_adjust < 0.05), ]

up_lnc_mRNA_090 <- lnc_mRNA_spearman_090[lnc_mRNA_spearman_090$lncRNA %in% up_lnc$gene, ]
down_lnc_mRNA_090 <- lnc_mRNA_spearman_090[lnc_mRNA_spearman_090$lncRNA %in% down_lnc$gene, ]
length(unique(up_lnc_mRNA_090$lncRNA))
length(unique(up_lnc_mRNA_090$mRNA))
length(unique(down_lnc_mRNA_090$lncRNA))
length(unique(down_lnc_mRNA_090$mRNA))

coexp_en_up <- enrichGO(gene = unique(up_lnc_mRNA_090$mRNA), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
unique(lnc_mRNA_spearman_090$mRNA)
coexp_up_res <- coexp_en_up@result
coexp_up_filt <- coexp_up_res[coexp_up_res$p.adjust < 0.05, ]
table(up_lnc_mRNA_090$lncRNA)

lnc_mRNA_spearman_085 <- lnc_mRNA_spearman[(abs(lnc_mRNA_spearman$correlation) >= 0.85 & lnc_mRNA_spearman$p_adjust < 0.05), ]

up_lnc_mRNA_085 <- lnc_mRNA_spearman_085[lnc_mRNA_spearman_085$lncRNA %in% up_lnc$gene, ]
up_lnc_mRNA_085 <- up_lnc_mRNA_085[order(up_lnc_mRNA_085$correlation), ]
rownames(up_lnc_mRNA_085) <- 1:nrow(up_lnc_mRNA_085)
up_lnc_mRNA_085 <- up_lnc_mRNA_085[-1, ]

length(unique(up_lnc_mRNA_085$lncRNA))
length(unique(up_lnc_mRNA_085$mRNA))

down_lnc_mRNA_085 <- lnc_mRNA_spearman_085[lnc_mRNA_spearman_085$lncRNA %in% down_lnc$gene, ]
length(unique(up_lnc_mRNA_085$lncRNA))
length(unique(up_lnc_mRNA_085$mRNA))
length(unique(down_lnc_mRNA_085$lncRNA))
length(unique(down_lnc_mRNA_085$mRNA))

length(unique(lnc_mRNA_spearman_085$mRNA))
length(unique(up_lnc_mRNA_085$mRNA))
length(unique(down_lnc_mRNA_085$mRNA))

coexp_en <- enrichGO(gene = unique(lnc_mRNA_spearman_085$mRNA), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
unique(lnc_mRNA_spearman_085$mRNA)
coexp_res <- coexp_en@result
coexp_filt <- coexp_res[coexp_res$p.adjust < 0.05, ]

coexp_en_up <- enrichGO(gene = unique(up_lnc_mRNA_085$mRNA), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
unique(up_lnc_mRNA_085$mRNA)
coexp_up_res <- coexp_en_up@result

coexp_en_down <- enrichGO(gene = unique(down_lnc_mRNA_085$mRNA), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
unique(down_lnc_mRNA_085$mRNA)
unique(lnc_mRNA_spearman_085$mRNA)
coexp_down_res <- coexp_en_down@result

up_terms <- c('defense response to virus', 'negative regulation of lymphocyte mediated immunity', 'positive regulation of interleukin-6 production', 'myeloid leukocyte activation')
down_terms <- c('mRNA splicing, via spliceosome', 'regulation of RNA splicing')

coexp_up_sel <- coexp_up_res[coexp_up_res$Description %in% up_terms, ]
coexp_down_sel <- coexp_down_res[coexp_down_res$Description %in% down_terms, ]

coexp_sel <- rbind(coexp_up_sel, coexp_down_sel)

coexp_sel$minus_log10_pvalue <- -log10(coexp_sel$pvalue)

coexp_sel$Description <- factor(coexp_sel$Description, levels = rev(c('defense response to virus', 'negative regulation of lymphocyte mediated immunity', 'positive regulation of interleukin-6 production', 'myeloid leukocyte activation', 'regulation of RNA splicing', 'mRNA splicing, via spliceosome')))

pdf('AgeAssociatedLncRNA_coexpmRNA_GOresult_barplot.pdf', height = 4, width = 8)
ggplot(coexp_sel, aes(x = minus_log10_pvalue, y = Description)) + geom_bar(position = "dodge", stat = "identity", width = 0.8)
dev.off()

sel_desc <- c('defense response to virus', 'negative regulation of lymphocyte mediated immunity', 'positive regulation of interleukin-6 production', 'myeloid leukocyte activation')

coexp_up_filt <- coexp_up_res[coexp_up_res$Description %in% sel_desc, ]

coexp_up_filt <- coexp_up_filt[, c(2, 8)]

coexp_up_split <- strsplit(coexp_up_filt$geneID, '/')
names(coexp_up_split) <- coexp_up_filt$Description
NAME <- names(coexp_up_split)
df <- data.frame(stringsAsFactors = F)
for (i in 1:4) {
  term_df <- data.frame(TERM = NAME[i], Gene = coexp_up_split[[i]])
  df <- rbind(term_df, df)
}
df$enrich_group <- 'GO_BP'
go_gene_df <- df

merge1 <- merge(up_lnc_mRNA_085, go_gene_df, by.x = 'mRNA', by.y = 'Gene', all.x = T)
merge1 <- merge1[!is.na(merge1$TERM), ]
merge1 <- unique(merge1)
unique(merge1$lncRNA)
merge1$lncRNA_direction <- 'upregulated'

go_gene_num <- coexp_up_res[, c(2, 9)]
colnames(go_gene_num)[2] <- 'Genenumber_of_each_GOterm'
merge1 <- merge(merge1, go_gene_num, by.x = 'TERM', by.y = 'Description', all.x = T)
length(unique(merge1$mRNA))
length(unique(merge1$lncRNA))

cancer_hallmark <- read.csv('./LncSEA/Cancer_Hallmark.csv', header = F, stringsAsFactors = F)
cancer_hallmark <- cancer_hallmark[, c(-1, -2)]
cancer_hallmark[cancer_hallmark == ''] <- NA
cancer_hallmark <- as.matrix(cancer_hallmark)
dim(cancer_hallmark) <- c(nrow(cancer_hallmark) * ncol(cancer_hallmark), 1)
cancer_hallmark <- as.data.frame(cancer_hallmark)
cancer_hallmark <- cancer_hallmark[!is.na(cancer_hallmark$V1), ]
cancer_hallmark <- as.data.frame(cancer_hallmark)

cancer_hallmark$Group <- 'Cancer_Hallmark'
colnames(cancer_hallmark) <- c('lncRNA_ID', 'Group')
cancer_hallmark <- unique(cancer_hallmark)

disease <- read.csv('./LncSEA/Disease.csv', header = F, stringsAsFactors = F)

disease <- disease[, c(-1, -2)]
disease[disease == ''] <- NA
disease <- as.matrix(disease)
dim(disease) <- c(nrow(disease) * ncol(disease), 1)
disease <- as.data.frame(disease)
disease <- disease[!is.na(disease$V1), ]
disease <- as.data.frame(disease)

disease$Group <- 'Disease'
colnames(disease) <- c('lncRNA_ID', 'Group')
disease <- unique(disease)

methylation <- read.csv('./LncSEA/Methylation_Pattern.csv', header = F, stringsAsFactors = F)

methylation <- methylation[, c(-1, -2)]
methylation[methylation == ''] <- NA
methylation <- as.matrix(methylation)
dim(methylation) <- c(nrow(methylation) * ncol(methylation), 1)
methylation <- as.data.frame(methylation)
methylation <- methylation[!is.na(methylation$V1), ]
methylation <- as.data.frame(methylation)

methylation$Group <- 'Methylation_Pattern'
colnames(methylation) <- c('lncRNA_ID', 'Group')
methylation <- unique(methylation)

microrna <- read.csv('./LncSEA/MicroRNA.csv', header = F, stringsAsFactors = F)

microrna <- microrna[, c(-1, -2)]
microrna[microrna == ''] <- NA
microrna <- as.matrix(microrna)
dim(microrna) <- c(nrow(microrna) * ncol(microrna), 1)
microrna <- as.data.frame(microrna)
microrna <- microrna[!is.na(microrna$V1), ]
microrna <- as.data.frame(microrna)

microrna$Group <- 'MicroRNA'
colnames(microrna) <- c('lncRNA_ID', 'Group')
microrna <- unique(microrna)

survival <- read.csv('./LncSEA/Survival.csv', header = F, stringsAsFactors = F)

survival <- survival[, c(-1, -2)]
survival[survival == ''] <- NA
survival <- as.matrix(survival)
dim(survival) <- c(nrow(survival) * ncol(survival), 1)
survival <- as.data.frame(survival)
survival <- survival[!is.na(survival$V1), ]
survival <- as.data.frame(survival)

survival$Group <- 'Survival'
colnames(survival) <- c('lncRNA_ID', 'Group')
survival <- unique(survival)

func_lnc_list <- rbind(cancer_hallmark, disease, methylation, microrna, survival)
func_lnc_list <- func_lnc_list[, c(2, 1)]
background_lnc_list <- lnc_list
background_lnc_list <- unique(background_lnc_list)
background_as <- background_lnc_list[background_lnc_list$LNCipedia_lncRNAlist %in% age_lnc$gene, ]
background_as <- data.frame(lncRNA_ID = background_as)
background_as$Group <- 'NA'
background_as <- background_as[, c(2, 1)]

background_nas <- background_lnc_list[!(background_lnc_list$LNCipedia_lncRNAlist %in% age_lnc$gene), ]
background_nas <- data.frame(lncRNA_ID = background_nas)
background_nas$Group <- 'NA'
background_nas <- background_nas[, c(2, 1)]

background_nas_sub <- background_nas[sample(1:nrow(background_nas), 34800), ]

func_lnc_with_bg <- rbind(func_lnc_list, background_nas_sub, background_as)

as_lnc_en <- enricher(gene = age_lnc$gene, TERM2GENE = func_lnc_with_bg, minGSSize = 0, maxGSSize = 30000)
as_lnc_res <- as_lnc_en@result

non_age_lnc <- read.table('./non_age_associated_lncRNA_selected2000_forLncSEA.txt', header = T, sep = '\t', stringsAsFactors = F)

nas_lnc_en <- enricher(gene = non_age_lnc$non_age_associated_lncRNA_selected, TERM2GENE = func_lnc_with_bg, minGSSize = 0, maxGSSize = 40000)
nas_lnc_res <- nas_lnc_en@result

enhancer <- read.csv('./LncSEA/Enhancer.csv', header = F, stringsAsFactors = F)
enhancer <- enhancer[, c(-1, -2)]
enhancer[enhancer == ''] <- NA
enhancer <- as.matrix(enhancer)
dim(enhancer) <- c(nrow(enhancer) * ncol(enhancer), 1)
enhancer <- as.data.frame(enhancer)
enhancer <- enhancer[!is.na(enhancer$V1), ]
enhancer <- as.data.frame(enhancer)

enhancer$Group <- 'Enhancer'
colnames(enhancer) <- c('lncRNA_ID', 'Group')
enhancer <- unique(enhancer)

super_enhancer <- read.csv('./LncSEA/Super_Enhancer.csv', header = F, stringsAsFactors = F)

super_enhancer <- super_enhancer[, c(-1, -2)]
super_enhancer[super_enhancer == ''] <- NA
super_enhancer <- as.matrix(super_enhancer)
dim(super_enhancer) <- c(nrow(super_enhancer) * ncol(super_enhancer), 1)
super_enhancer <- as.data.frame(super_enhancer)
super_enhancer <- super_enhancer[!is.na(super_enhancer$V1), ]
super_enhancer <- as.data.frame(super_enhancer)

super_enhancer$Group <- 'Super_Enhancer'
colnames(super_enhancer) <- c('lncRNA_ID', 'Group')
super_enhancer <- unique(super_enhancer)

tf <- read.csv('./LncSEA/Transcription_Factor.csv', header = F, stringsAsFactors = F)

tf <- tf[, c(-1, -2)]
tf[tf == ''] <- NA
tf <- as.matrix(tf)
dim(tf) <- c(nrow(tf) * ncol(tf), 1)
tf <- as.data.frame(tf)
tf <- tf[!is.na(tf$V1), ]
tf <- as.data.frame(tf)

tf$Group <- 'Transcription_Factor'
colnames(tf) <- c('lncRNA_ID', 'Group')
tf <- unique(tf)

rbp <- read.csv('./LncSEA/RNA_Binding_Protein.csv', header = F, stringsAsFactors = F)

rbp <- rbp[, c(-1, -2)]
rbp[rbp == ''] <- NA
rbp <- as.matrix(rbp)
dim(rbp) <- c(nrow(rbp) * ncol(rbp), 1)
rbp <- as.data.frame(rbp)
rbp <- rbp[!is.na(rbp$V1), ]
rbp <- as.data.frame(rbp)

rbp$Group <- 'RNA_Binding_Protein'
colnames(rbp) <- c('lncRNA_ID', 'Group')
rbp <- unique(rbp)

chromatin <- read.csv('./LncSEA/Accessible_Chromatin.csv', header = F, stringsAsFactors = F)

chromatin <- chromatin[, c(-1, -2)]
chromatin[chromatin == ''] <- NA
chromatin <- as.matrix(chromatin)
dim(chromatin) <- c(nrow(chromatin) * ncol(chromatin), 1)
chromatin <- as.data.frame(chromatin)
chromatin <- chromatin[!is.na(chromatin$V1), ]
chromatin <- as.data.frame(chromatin)

chromatin$Group <- 'Accessible_Chromatin'
colnames(chromatin) <- c('lncRNA_ID', 'Group')
chromatin <- unique(chromatin)

sum(tf$lncRNA_ID %in% enhancer$lncRNA_ID)
sum(chromatin$lncRNA_ID %in% enhancer$lncRNA_ID)

func_lnc_list <- rbind(enhancer, super_enhancer, tf, rbp, chromatin)
func_lnc_list <- func_lnc_list[, c(2, 1)]

background_lnc_list <- lnc_list
background_lnc_list <- unique(background_lnc_list)
background_as <- background_lnc_list[background_lnc_list$LNCipedia_lncRNAlist %in% age_lnc$gene, ]
background_as <- data.frame(lncRNA_ID = background_as)
background_as$Group <- 'NA'
background_as <- background_as[, c(2, 1)]

background_nas <- background_lnc_list[!(background_lnc_list$LNCipedia_lncRNAlist %in% age_lnc$gene), ]
background_nas <- data.frame(lncRNA_ID = background_nas)
background_nas$Group <- 'NA'
background_nas <- background_nas[, c(2, 1)]

background_nas_sub <- background_nas[sample(1:nrow(background_nas), 0), ]

func_lnc_with_bg <- rbind(func_lnc_list, background_nas_sub, background_as)

as_lnc_en <- enricher(gene = age_lnc$gene, TERM2GENE = func_lnc_with_bg, minGSSize = 0, maxGSSize = 30000)
as_lnc_res <- as_lnc_en@result

non_age_lnc <- read.table('./non_age_associated_lncRNA_selected2000_forLncSEA.txt', header = T, sep = '\t', stringsAsFactors = F)

nas_lnc_en <- enricher(gene = non_age_lnc$non_age_associated_lncRNA_selected, TERM2GENE = func_lnc_with_bg, minGSSize = 0, maxGSSize = 30000)
nas_lnc_res <- nas_lnc_en@result

upreg_res <- read.table('./age_associated_lncRNA_LncSEA_5UpstreamRegulation_geneSet_EnrichResult.txt', header = T, sep = ' ', stringsAsFactors = F)
func_res <- read.table('./age_associated_lncRNA_LncSEA_5functional_geneSet_EnrichResult.txt', header = T, sep = '\t', stringsAsFactors = F)

all_res <- rbind(upreg_res, func_res)
all_res <- all_res[, c(2, 6, 9)]
all_res$group <- c(rep('Regulational', 5), rep('Functional', 5))
all_res$minus_log10 <- -log10(all_res$p.adjust)
all_res$Description <- factor(all_res$Description, levels = all_res$Description)
all_res$minus_log10[all_res$minus_log10 >= 70] <- 70

pdf('LncSEA_EnrichmentAnalysis_BubblePlot.pdf', height = 4, width = 7)
ggplot(all_res, aes(group, Description)) + geom_point(aes(color = minus_log10, size = Count)) + labs(x = NULL, y = NULL) + scale_color_gradient(low = "#A0DA9A", high = "#02471D", limit = c(1.30103, 70))
dev.off()

length(colnames(tpm_nochi))

subcell <- read.csv('./LncSEA/Subcellular_Location.csv', header = F, stringsAsFactors = F)

subcell <- subcell[1:10, ]

nucleus <- subcell[1, ]
nucleus <- as.character(nucleus)
nucleus <- nucleus[c(-1, -2)]
nucleus <- nucleus[nucleus != '']
nucleus <- data.frame(lncRNA = nucleus)
nucleus$group <- 'Nucleus'

cytosol <- subcell[2, ]
cytosol <- as.character(cytosol)
cytosol <- cytosol[c(-1, -2)]
cytosol <- cytosol[cytosol != '']
cytosol <- data.frame(lncRNA = cytosol)
cytosol$group <- 'Cytosol'

cytoplasm <- subcell[3, ]
cytoplasm <- as.character(cytoplasm)
cytoplasm <- cytoplasm[c(-1, -2)]
cytoplasm <- cytoplasm[cytoplasm != '']
cytoplasm <- data.frame(lncRNA = cytoplasm)
cytoplasm$group <- 'Cytoplasm'

exosome <- subcell[4, ]
exosome <- as.character(exosome)
exosome <- exosome[c(-1, -2)]
exosome <- exosome[exosome != '']
exosome <- data.frame(lncRNA = exosome)
exosome$group <- 'Exosome'

ribosome <- subcell[5, ]
ribosome <- as.character(ribosome)
ribosome <- ribosome[c(-1, -2)]
ribosome <- ribosome[ribosome != '']
ribosome <- data.frame(lncRNA = ribosome)
ribosome$group <- 'Ribosome'

er <- subcell[6, ]
er <- as.character(er)
er <- er[c(-1, -2)]
er <- er[er != '']
er <- data.frame(lncRNA = er)
er$group <- 'Endoplasmic_reticulum'

mito <- subcell[7, ]
mito <- as.character(mito)
mito <- mito[c(-1, -2)]
mito <- mito[mito != '']
mito <- data.frame(lncRNA = mito)
mito$group <- 'Mitochondrion'

circ <- subcell[8, ]
circ <- as.character(circ)
circ <- circ[c(-1, -2)]
circ <- circ[circ != '']
circ <- data.frame(lncRNA = circ)
circ$group <- 'Circulating'

nucleolus <- subcell[9, ]
nucleolus <- as.character(nucleolus)
nucleolus <- nucleolus[c(-1, -2)]
nucleolus <- nucleolus[nucleolus != '']
nucleolus <- data.frame(lncRNA = nucleolus)
nucleolus$group <- 'Nucleolus'

nucleoplasm <- subcell[10, ]
nucleoplasm <- as.character(nucleoplasm)
nucleoplasm <- nucleoplasm[c(-1, -2)]
nucleoplasm <- nucleoplasm[nucleoplasm != '']
nucleoplasm <- data.frame(lncRNA = nucleoplasm)
nucleoplasm$group <- 'Nucleoplasm'

subcell_info <- rbind(nucleus, cytosol, cytoplasm, exosome, ribosome, er, mito, circ, nucleolus, nucleoplasm)
table(subcell_info$group)
length(unique(subcell_info$lncRNA))

as_subcell <- subcell_info[subcell_info$lncRNA %in% age_lnc$gene, ]

subcell <- read.csv('./LncSEA/Subcellular_Location.csv', header = F, stringsAsFactors = F)

cytoplasm_iloc <- subcell[11, ]
cytoplasm_iloc <- as.character(cytoplasm_iloc)
cytoplasm_iloc <- cytoplasm_iloc[c(-1, -2)]
cytoplasm_iloc <- cytoplasm_iloc[cytoplasm_iloc != '']
cytoplasm_iloc <- data.frame(lncRNA = cytoplasm_iloc)
cytoplasm_iloc$group <- 'Cytoplasm_iLoc'

exosome_iloc <- subcell[12, ]
exosome_iloc <- as.character(exosome_iloc)
exosome_iloc <- exosome_iloc[c(-1, -2)]
exosome_iloc <- exosome_iloc[exosome_iloc != '']
exosome_iloc <- data.frame(lncRNA = exosome_iloc)
exosome_iloc$group <- 'Exosome_iLoc'

nucleus_iloc <- subcell[13, ]
nucleus_iloc <- as.character(nucleus_iloc)
nucleus_iloc <- nucleus_iloc[c(-1, -2)]
nucleus_iloc <- nucleus_iloc[nucleus_iloc != '']
nucleus_iloc <- data.frame(lncRNA = nucleus_iloc)
nucleus_iloc$group <- 'Nucleus_iLoc'

ribosome_iloc <- subcell[14, ]
ribosome_iloc <- as.character(ribosome_iloc)
ribosome_iloc <- ribosome_iloc[c(-1, -2)]
ribosome_iloc <- ribosome_iloc[ribosome_iloc != '']
ribosome_iloc <- data.frame(lncRNA = ribosome_iloc)
ribosome_iloc$group <- 'Ribosome_iLoc'

iloc <- rbind(cytoplasm_iloc, exosome_iloc, nucleus_iloc, ribosome_iloc)

as_iloc <- iloc[iloc$lncRNA %in% age_lnc$gene, ]
length(unique(as_iloc$lncRNA))
table(as_iloc$group)

as_iloc_bar <- as.data.frame(table(as_iloc$group))
as_iloc_bar <- as_iloc_bar[order(as_iloc_bar$Freq), ]

as_iloc_bar$Var1 <- factor(as_iloc_bar$Var1, levels = rev(as_iloc_bar$Var1))

pdf('age_associated_lncRNA_subcellularLocation_byiLocFromLncSEA_barplot.pdf', height = 4, width = 6)
ggplot(as_iloc_bar, aes(x = Var1, y = Freq)) + geom_bar(position = 'dodge', stat = "identity", width = 0.4)
dev.off()