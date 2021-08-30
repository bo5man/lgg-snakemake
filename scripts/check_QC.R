###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(openxlsx))
    suppressMessages(library(minfi))
    suppressMessages(library(RSpectra))
    suppressMessages(library(dplyr))
} else{
    library(openxlsx)
    library(minfi)
    library(RSpectra)
    library(dplyr)
}

##################
# input file paths
# sentrixs_full_cohort <- snakemake@input[['sentrix']]
# pRGset <- snakemake@input[["RGset"]]
pMset_full <- snakemake@input[["Mset_full"]]
pBetas_full <- snakemake@input[["betas_full"]]
pDFclinical_full_cohort <- snakemake@input[['DFclinical_full_cohort']]
pDFclinical_full_gliomas <- snakemake@input[['DFclinical_full_gliomas']]
pDFclinical_full_inhouse <- snakemake@input[['DFclinical_full_inhouse']]

# output file paths
pqcMset_full_cohort <- snakemake@output[["qcMset_full_cohort"]]
pqcplotMset_full_cohort <- snakemake@output[["qcplotMset_full_cohort"]]
pdensityplotMset_full_cohort <- snakemake@output[['densityplotMset_full_cohort']]
pDFclinical_cohort <- snakemake@output[['DFclinical_cohort']]
pDFclinical_gliomas <- snakemake@output[['DFclinical_gliomas']]
pDFclinical_inhouse <- snakemake@output[['DFclinical_inhouse']]
pDFclinical_cohort_csv <- snakemake@output[['DFclinical_cohort_csv']]
pDFclinical_inhouse_csv <- snakemake@output[['DFclinical_inhouse_csv']]
pBetas_cohort<- snakemake@output[["betas_cohort"]]
pBetas_gliomas <- snakemake@output[["betas_gliomas"]]
pBetas_inhouse <- snakemake@output[["betas_inhouse"]]

# parameters
# # parameters
# pinhouse <- './../../LGG_Methylation/Yongsoo VUMC methylation arrat samples sep 2020.xlsx'
# pcohort <- './../../LGG_Methylation/Overview_oligodendroglioma_samples_long_vs_short_PFS.xlsx'
pcohort <- snakemake@params[['overview']]
dir_full_cohort <-  snakemake@params[['dir_full_cohort']]
dir_qc <- snakemake@params[["dir_qc"]]
if(!dir.exists(dir_qc))(dir.create(dir_qc))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# # # read input
# RGset <- readRDS(pRGset)
Mset_full <- readRDS(pMset_full)
betas_full <- readRDS(pBetas_full)
DFclinical_full_cohort <- readRDS(pDFclinical_full_cohort)
DFclinical_full_gliomas <- readRDS(pDFclinical_full_gliomas)
DFclinical_full_inhouse <- readRDS(pDFclinical_full_inhouse)

Mset_full_cohort <- Mset_full[,rownames(DFclinical_full_cohort)]
betas_full_cohort <- betas_full[rownames(DFclinical_full_inhouse),]


# # # QC for Mset
qcMset_full_cohort <- minfiQC(Mset_full_cohort, verbose=T)
message('saving quality control of Mset of cohort samples in ', pqcMset_full_cohort)
saveRDS(qcMset_full_cohort, file = pqcMset_full_cohort)

qc <- qcMset_full_cohort$qc
#plotQC
#?plotQC
badSampleCutoff <- 12
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < badSampleCutoff))
sentrixBadQCCohort <- rownames(qc)[whichBad]
namesBadQCCohort <- DFclinical_full_cohort[sentrixBadQCCohort,'ID']

message('saving plot of quality control of Mset of full cohort samples in ', pqcplotMset_full_cohort)
png(file = pqcplotMset_full_cohort, width = 1440, height = 1280) #width = 1846,height = 991
plot(qc$mMed, qc$uMed, xlim = c(9, 15), ylim = c(9, 15),
    xaxt = "n", yaxt = "n", xlab = "Meth median intensity (log2)",
    ylab = "Unmeth median intensity (log2)", col = ifelse(1:nrow(qc) %in%
        whichBad, "red", "black"))
axis(side = 1, at = c(10, 12, 14))
axis(side = 2, at = c(10, 12, 14))
abline(badSampleCutoff * 2, -1, lty = 2)
text(11,13 + 0.35, 'uMed = 24 - mMed')
title('Quality control on median intensities of Methylation set for samples of the full cohort')
text(qc$mMed, qc$uMed + 0.25, labels = DFclinical_full_cohort[rownames(qc),'LSS'], col = ifelse(1:nrow(qc) %in% whichBad, 'red', 'black'))
dev.off()

message('Bad quality samples found by minfiQC: ', paste0(namesBadQCCohort,' '))
# Bad quality samples found: LSS16  LSS19  LSS11  LSS27

# # # densityplot for Mset of full cohort
message('saving densityplot of Mset of full cohort samples in ', pdensityplotMset_full_cohort)
png(file = pdensityplotMset_full_cohort, width = 1440, height = 1280) #width = 1846,height = 991
par(mfrow=c(2,1))
densityPlot(Mset_full_cohort, sampGroups=paste(DFclinical_full_cohort[colnames(Mset_full_cohort),'TypeSurvival'], 'survivor'))
title('Densityplot of Mset of samples of the full cohort, categorized by survival')
densityPlot(Mset_full_cohort, sampGroups=ifelse(1:ncol(Mset_full_cohort) %in% whichBad, 'bad','good'))
title('Densityplot of Mset of samples of the full cohort, categorized by quality')
dev.off()



# # # read cohort Overview_oligodendroglioma_samples_long_vs_short_PFS

cohort_full <- read.xlsx(pcohort, sep = '_')
cohort_full <- cohort_full[!is.na(cohort_full$SENTRIX_ID),]
sentrixBadQCCohort <-  sentrixBadQCCohort
# # # samples with an q37/total ratio less than 0.25 will not be used
badSampleCutoffQDNAseq <- 0.25
thresQDNAseq <- cohort_full$QDNAseq_Ratio_q37Total < badSampleCutoffQDNAseq
sentrixBadQDNASeqRatio <- cohort_full$SENTRIX_ID[thresQDNAseq]
# # # samples that did not have a 1p/19q co-deletion in the methylation cnv will not be used
typesBad1p19qCohort <- grep('1p/19q positive', cohort_full$Methadone_profile_Own_Script_Result_comments, value=T, invert = T)
sentrixBad1p19qCohort <- cohort_full$SENTRIX_ID[cohort_full$Methadone_profile_Own_Script_Result_comments %in% typesBad1p19qCohort]
namesBad1p19qCohort <- cohort_full$`Sample_ID_Long_vs_Short_Survivor_(LSS)`[cohort_full$Methadone_profile_Own_Script_Result_comments %in% typesBad1p19qCohort]

message('Samples without 1p/19q found: ', paste0(namesBad1p19qCohort,' '))
# Samples without 1p/19q found: LSS11 LSS16 LSS18 LSS19 LSS21 LSS23 LSS24

sentrixBadCohort <- union(sentrixBadQCCohort,sentrixBad1p19qCohort)
namesBadCohort <- union(namesBadQCCohort,namesBad1p19qCohort)

message('Name IDs of removed samples for cohort: ', paste0(sort(namesBadCohort),' '))
# Name IDs of removed samples for cohort: LSS11 LSS16 LSS18 LSS19 LSS21 LSS23 LSS24 LSS27

# # # filter full clinical data frame to define clinical data frame

DFclinical_cohort <- DFclinical_full_cohort[!rownames(DFclinical_full_cohort) %in% sentrixBadCohort,]
message('saving clinical dataframe for good samples of 1p/19q of cohort in ', pDFclinical_cohort)
saveRDS(DFclinical_cohort, file = pDFclinical_cohort)

message('saving clinical dataframe as cnv of good samples of 1p/19q cohort in ', pDFclinical_cohort_csv)
write.csv2(DFclinical_cohort, file = pDFclinical_cohort_csv)


# # # do the same for DFclinical_gliomas
DFclinical_gliomas <- DFclinical_full_gliomas[!rownames(DFclinical_full_gliomas) %in% sentrixBadCohort,]
message('saving clinical dataframe for good samples of 1p/19q of cohort and gliomas in ', pDFclinical_gliomas)
saveRDS(DFclinical_gliomas, file = pDFclinical_gliomas)

# # # do the same for DFclinical_inhouse
DFclinical_inhouse <- DFclinical_full_inhouse[!rownames(DFclinical_full_inhouse) %in% sentrixBadCohort,]
message('saving clinical dataframe for good samples of 1p/19q of cohort and inhouse in ', pDFclinical_inhouse)
saveRDS(DFclinical_inhouse, file = pDFclinical_inhouse)

message('saving clinical dataframe as cnv of good samples of 1p/19q cohort and inhouse in ', pDFclinical_inhouse_csv)
write.csv2(DFclinical_inhouse, file = pDFclinical_inhouse_csv)

# # define betas for cohort

i <- intersect(rownames(DFclinical_cohort),rownames(betas_full))
betas_cohort <- betas_full[i,]

message('saving betas for good samples of cohort in ', pBetas_cohort)
saveRDS(betas_cohort, file = pBetas_cohort)


# # define betas for cohort and gliomas

i <- intersect(rownames(DFclinical_gliomas),rownames(betas_full))
betas_gliomas <- betas_full[i,]

message('saving betas for good samples of cohort and gliomas in ', pBetas_gliomas)
saveRDS(betas_gliomas, file = pBetas_gliomas)

# # define betas for cohort and inhouse

i <- intersect(rownames(DFclinical_inhouse),rownames(betas_full))
betas_inhouse <- betas_full[i,]

message('saving betas for good samples of cohort and inhouse in ', pBetas_inhouse)
saveRDS(betas_inhouse, file = pBetas_inhouse)

sink()
