# Christine Hwang, Evan Leung
# FCBB: Final Project
# Cancer Immunotherapy (Topic 7)

library(BSgenome.Hsapiens.1000genomes.hs37d5)
library('SomaticSignatures')
library('SomaticCancerAlterations')
library('TCGAbiolinks')

## preprocessing TCGA WES studies
sca_data = unlist(scaLoadDatasets())

sca_data$study = factor(gsub("(.*)_(.*)", "\\1", toupper(names(sca_data))))
sca_data = unname(subset(sca_data, Variant_Type %in% "SNP"))
sca_data = keepSeqlevels(sca_data, hsAutosomes(), pruning.mode = "coarse")
sca_data = sca_data[sca_data$study %in% c('HNSC'), ]

sca_vr = VRanges(
    seqnames = seqnames(sca_data),
    ranges = ranges(sca_data),
    ref = sca_data$Reference_Allele,
    alt = sca_data$Tumor_Seq_Allele2,
    sampleNames = sca_data$Patient_ID,
    seqinfo = seqinfo(sca_data),
    study = sca_data$study)

## extracting motifs
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
sca_mm = motifMatrix(sca_motifs, normalize = TRUE)
sca_mm_new <- sca_mm[ , colSums(is.na(sca_mm)) < nrow(sca_mm)]

## inferring somatic signatures
# n_sigs = 2:8
# gof_nmf = assessNumberSignatures(sca_mm_new, n_sigs, nReplicates = 5)
sigs_nmf = identifySignatures(sca_mm_new, 8, pcaDecomposition)
# sigs_nmf

## making association data
association_data = samples(sigs_nmf)
HNSC_infiltration <- read.delim('/Users/chrhwang/Documents/School/Spring22/FCBB/final_project/head_and_neck_squamous_cell_carcinoma_RNAseqV2.txt')
association_df <- data.frame(
	Patient_ID = names(association_data[, 'S1']),
	S1 = unname(association_data[, 'S1']),
	S2 = unname(association_data[, 'S2']),
	S3 = unname(association_data[, 'S3']),
	S4 = unname(association_data[, 'S4']),
	S5 = unname(association_data[, 'S5']),
	S6 = unname(association_data[, 'S6']),
	S7 = unname(association_data[, 'S7']),
	S8 = unname(association_data[, 'S8']),
	Infiltration_score =  rep(NA, length(names(association_data[, 'S1']))),
	Infiltration_l_h = rep('None', length(names(association_data[, 'S1']))))
for (i in HNSC_infiltration$ID) {
	if (substring(i, 1, 12) %in% association_df$Patient_ID) {
		association_df[association_df$Patient_ID == substring(i, 1, 12), ]$Infiltration_score = HNSC_infiltration[HNSC_infiltration$ID == i, ]$Immune_score
	}
}
	
## checking for association
association_df <- na.omit(association_df)
med = median(association_df$Infiltration_score)
for (i in association_df$Infiltration_score) {
	if (i <= med) {
		association_df[association_df$Infiltration_score == i, ]$Infiltration_l_h = 'Low Infiltration Score'
	}
	else {
		association_df[association_df$Infiltration_score == i, ]$Infiltration_l_h = 'High Infiltration Score'
	}
}

chisq <- chisq.test(association_df)