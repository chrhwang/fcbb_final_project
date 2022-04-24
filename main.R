# Christine Hwang, Evan Leung
# FCBB: Final Project
# Cancer Immunotherapy (Topic 7)

library(BSgenome.Hsapiens.1000genomes.hs37d5)
library('GenomicRanges')
library('SomaticSignatures')
library('SomaticCancerAlterations')
library('TCGAbiolinks')

# preprocessing TCGA WES studies
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

# extracting motifs
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
sca_mm = motifMatrix(sca_motifs, normalize = TRUE)
sca_mm_new <- sca_mm[ , colSums(is.na(sca_mm)) < nrow(sca_mm)]

# inferring somatic signatures
# n_sigs = 2:8
# gof_nmf = assessNumberSignatures(sca_mm_new, n_sigs, nReplicates = 5)
sigs_nmf = identifySignatures(sca_mm_new, 8, pcaDecomposition)
# sigs_nmf
samples(sigs_nmf)