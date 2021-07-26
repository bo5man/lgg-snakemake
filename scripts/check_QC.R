###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
    suppressMessages(library(RSpectra))
    suppressMessages(library(dplyr))
} else{
    library(minfi)
    library(RSpectra)
    library(dplyr)
}

##################
# input file paths
sentrixs_full_cohort <- snakemake@input[['sentrix']]
pRGset <- snakemake@input[["RGset"]]
pMset <- snakemake@input[["Mset"]]
pBetas <- snakemake@input[["betas"]]
pDFclinical_full_cohort <- snakemake@input[['DFclinical_full_cohort']]
pDFclinical_full_gliomas <- snakemake@input[['DFclinical_full_gliomas']]
pDFclinical_full_inhouse <- snakemake@input[['DFclinical_full_inhouse']]
# sentrixs_full_cohort = list.files(path='./../../LGG_Methylation/data/utrecht_methy/', pattern='*_Grn.idat')

# output file paths
pqcMset_full_cohort <- snakemake@output[["qcMset_full_cohort"]]
pqcplotMset_full_cohort <- snakemake@output[["qcplotMset_full_cohort"]]
pdensityplotMset_full_cohort <- snakemake@output[['densityplotMset_full_cohort']]
pDFclinical_cohort <- snakemake@output[['DFclinical_cohort']]
pDFclinical_gliomas <- snakemake@output[['DFclinical_gliomas']]
pDFclinical_inhouse <- snakemake@output[['DFclinical_inhouse']]
pBetas_cohort<- snakemake@output[["betas_cohort"]]
pBetas_gliomas <- snakemake@output[["betas_gliomas"]]
pBetas_inhouse <- snakemake@output[["betas_inhouse"]]

# parameters
dir_full_cohort <-  snakemake@params[['dir_full_cohort']]
dir_qc <- snakemake@params[["dir_qc"]]
if(!dir.exists(dir_qc))(dir.create(dir_qc))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

# # # read input
RGset <- readRDS(pRGset)
Mset <- readRDS(pMset)
betas <- readRDS(pBetas)
DFclinical_full_cohort <- readRDS(pDFclinical_full_cohort)
DFclinical_full_gliomas <- readRDS(pDFclinical_full_gliomas)
DFclinical_full_inhouse <- readRDS(pDFclinical_full_inhouse)

# # # filter RGset and Mset on full cohort samples
sentrixs_full_cohort <- gsub(dir_full_cohort, '', sentrixs_full_cohort)
sentrixs_full_cohort <- gsub('_Grn.idat', '', sentrixs_full_cohort)
RGset_full_cohort <- RGset[,sentrixs_full_cohort]
Mset_full_cohort <- Mset[,sentrixs_full_cohort]
TypeSurvival <- paste(DFclinical_full_cohort[rownames(DFclinical_full_cohort) %in% sentrixs_full_cohort,'TypeSurvival'], 'survivor')

# # # QC for Mset
qcMset_full_cohort <- minfiQC(Mset_full_cohort, verbose=T)
message('saving quality control of Mset of cohort samples in ', pqcMset_full_cohort)
saveRDS(qcMset_full_cohort, file = pqcMset_full_cohort)

qc = qcMset_full_cohort$qc
#plotQC
#?plotQC
badSampleCutoff <- 12
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < badSampleCutoff))
sentrixBadQCCohort <- rownames(qc)[whichBad]
namesBad <- DFclinical_full_cohort[sentrixBadQCCohort,'LSS']

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

message('Bad quality samples found: ', paste0(namesBad,' '))
# Bad quality samples found: LSS16  LSS19  LSS11  LSS27
# print(DFclinical_full_cohort[sentrixBad,])

# # # densityplot for Mset of full cohort
message('saving densityplot of Mset of full cohort samples in ', pdensityplotMset_full_cohort)
png(file = pdensityplotMset_full_cohort, width = 1440, height = 1280) #width = 1846,height = 991
par(mfrow=c(2,1))
densityPlot(Mset_full_cohort, sampGroups=TypeSurvival)
title('Densityplot of Mset of samples of the full cohort, categorized by survival')
densityPlot(Mset_full_cohort, sampGroups=ifelse(1:ncol(Mset_full_cohort) %in% whichBad, 'bad','good'))
title('Densityplot of Mset of samples of the full cohort, categorized by quality')
dev.off()


typesBad1p19qCohort <- grep('1p/19q', DFclinical_full_cohort$Type, value=T, invert = T)
sentrixBad1p19qCohort <- DFclinical_full_cohort$Sentrix_ID[DFclinical_full_cohort$Type %in% typesBad1p19qCohort]
sentrixBadCohort <- union(sentrixBadQCCohort,sentrixBad1p19qCohort)
sentrixGoodCohort <- setdiff(DFclinical_full_cohort$Sentrix_ID, sentrixBadCohort)
DFclinical_cohort <- DFclinical_full_cohort[sentrixGoodCohort,]

message('saving clinical dataframe for good samples of 1p/19q of cohort in ', pDFclinical_cohort)
saveRDS(DFclinical_cohort, file = pDFclinical_cohort)

sentrixGoodGliomas <- setdiff(DFclinical_full_gliomas$Sentrix_ID,sentrixBadCohort)
DFclinical_gliomas <- DFclinical_full_gliomas[sentrixGoodGliomas,]
message('saving clinical dataframe for good samples of 1p/19q of cohort and gliomas in ', pDFclinical_gliomas)
saveRDS(DFclinical_gliomas, file = pDFclinical_gliomas)

sentrixGoodInhouse <- setdiff(DFclinical_full_inhouse$Sentrix_ID,sentrixBadCohort)
DFclinical_inhouse <- DFclinical_full_inhouse[sentrixGoodInhouse,]
message('saving clinical dataframe for good samples of 1p/19q of cohort and inhouse in ', pDFclinical_inhouse)
saveRDS(DFclinical_inhouse, file = pDFclinical_inhouse)


# # define betas for cohort

i <- intersect(rownames(DFclinical_cohort),rownames(betas))
betas_cohort <- betas[i,]

message('saving betas for good samples of cohort in ', pBetas_cohort)
saveRDS(betas_cohort, file = pBetas_cohort)

# # define betas for cohort and gliomas

i <- intersect(rownames(DFclinical_gliomas),rownames(betas))
betas_gliomas <- betas[i,]

message('saving betas for good samples of cohort and gliomas in ', pBetas_gliomas)
saveRDS(betas_gliomas, file = pBetas_gliomas)

# # define betas for cohort and inhouse

i <- intersect(rownames(DFclinical_inhouse),rownames(betas))
betas_inhouse <- betas[i,]

message('saving betas for good samples of cohort and inhouse in ', pBetas_inhouse)
saveRDS(betas_inhouse, file = pBetas_inhouse)
