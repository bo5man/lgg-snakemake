###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(openxlsx))
    suppressMessages(library(dplyr))
    suppressMessages(library(stringr))
    suppressMessages(library(gplots))
    suppressMessages(library(minfi))
    suppressMessages(library(purrr))
    suppressMessages(library(RSpectra))
    suppressMessages(library(Rtsne))
    suppressMessages(library(ggplot2))
} else{
    library(openxlsx)
    library(dplyr)
    library(stringr)
    library(gplots)
    library(minfi)
    library(purrr)
    library(RSpectra)
    library(Rtsne)
    library(ggplot2)
}

##################
# # input file paths
pMset <- snakemake@input[["Mset"]]
pbetas <- snakemake@input[["betas"]]
# pbetas <- './../results/betas/betas.RDS'

# # output file paths
pDFclinical_full_cohort <-      snakemake@output[["DFclinical_full_cohort"]]
pDFclinical_full_gliomas <-      snakemake@output[["DFclinical_full_gliomas"]]
pDFclinical_full_inhouse <-      snakemake@output[["DFclinical_full_inhouse"]]
# pDFclinical_full_gliomas <-      './../results/DFclinical_full_gliomas.RDS'
pDFclinical_full_cohort_csv <-      snakemake@output[["DFclinical_full_cohort_csv"]]
pDFclinical_full_inhouse_csv <-      snakemake@output[["DFclinical_full_inhouse_csv"]]
pMset_full <- snakemake@output[["Mset_full"]]
pbetas_full <- snakemake@output[["betas_full"]]

# # parameters
# pinhouse <- './../../LGG_Methylation/Yongsoo VUMC methylation arrat samples sep 2020.xlsx'
# pcohort <- './../../LGG_Methylation/Overview_oligodendroglioma_samples_long_vs_short_PFS.xlsx'
pinhouse <- snakemake@params[['inhouse']]
pcohort <- snakemake@params[['overview']]

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# read input
Mset <- readRDS(pMset)
betas <- readRDS(pbetas)

# # # read clinical data files
inhouse_full <- read.xlsx(pinhouse, sep = "_")
cohort_full <- read.xlsx(pcohort, sep = '_')

# # # select inhouse data and full cohort data and create data frame sDF_inhouse
inhouse <- inhouse_full[,c('Sentrix_ID', 'Methylatie_array_uitkomst_schoon', 'Classifier_methylation_subclass_(only_for_predicted_methylation_class_families,_cut-off_0.5)_schoon')]
# # # select cohort data
# cohort <- cohort_full[,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Shallow_Sequencing_Result', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]
cohort <- cohort_full[!is.na(cohort_full$SENTRIX_ID),c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Methylation_classes','Methylation_subclasses', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]

Selected1p19q <- unique(str_subset(inhouse$`Classifier_methylation_subclass_(only_for_predicted_methylation_class_families,_cut-off_0.5)_schoon`, '1p/19q'))

iFullCohort <- seq_along(cohort$SENTRIX_ID)
iSelected1p19q <- which(inhouse$`Classifier_methylation_subclass_(only_for_predicted_methylation_class_families,_cut-off_0.5)_schoon` %in% Selected1p19q)
iOtherSelectedGliomas <- setdiff(seq_along(inhouse$`Classifier_methylation_subclass_(only_for_predicted_methylation_class_families,_cut-off_0.5)_schoon`), iSelected1p19q)


sDF_full_inhouse <- data.frame(
                         ID = c(cohort$`Sample_ID_Long_vs_Short_Survivor_(LSS)`[iFullCohort],
                                 paste0('1p19qI', seq_along(inhouse$Methylatie_array_uitkomst_schoon[iSelected1p19q])),
                                 paste0('OG', seq_along(inhouse$Methylatie_array_uitkomst_schoon[iOtherSelectedGliomas]))
                                 ),
                         Cohort_ID = c(rep('FullCohort',length(cohort$SENTRIX_ID[iFullCohort])),
                                 rep('Selected1p19q',length(inhouse$Sentrix_ID[iSelected1p19q])),
                                 rep('OtherSelectedGliomas',length(inhouse$Sentrix_ID[iOtherSelectedGliomas]))
                                 ),
                         Sentrix_ID = c(cohort$SENTRIX_ID[iFullCohort],
                                 inhouse$Sentrix_ID[iSelected1p19q],
                                 inhouse$Sentrix_ID[iOtherSelectedGliomas]
                                 ),
                         Type = c(ifelse(is.na(cohort$Methylation_subclasses[iFullCohort]),cohort$Methylation_classes[iFullCohort],paste(cohort$Methylation_classes[iFullCohort], cohort$Methylation_subclasses[iFullCohort],sep=', ')),
                                 ifelse(is.na(inhouse[iSelected1p19q,3]), inhouse[iSelected1p19q,2], paste(inhouse[iSelected1p19q,2],inhouse[iSelected1p19q,3], sep=', ')),
                                 ifelse(is.na(inhouse[iOtherSelectedGliomas,3]), inhouse[iOtherSelectedGliomas,2], paste(inhouse[iOtherSelectedGliomas,2],inhouse[iOtherSelectedGliomas,3], sep=', '))
                                 ),
                        #TypeSurvival = c(paste('cohort', cohort$`long_or_short_(>100_months_-_<40_months)`[iFullCohort], 'survivor'),
                         TypeSurvival = c(paste(cohort$`long_or_short_(>100_months_-_<40_months)`[iFullCohort], 'survivor'),
                                 rep('Selected1p19q',length(inhouse$Sentrix_ID[iSelected1p19q])),
                                 rep('OtherSelectedGliomas',length(inhouse$Sentrix_ID[iOtherSelectedGliomas]))
                                 ),
                         SexType = c(cohort$`Sex_(M/F)`[iFullCohort],
                                 rep('Selected1p19q',length(inhouse$Sentrix_ID[iSelected1p19q])),
                                 rep('OtherSelectedGliomas',length(inhouse$Sentrix_ID[iOtherSelectedGliomas]))
                                 )
)
rownames(sDF_full_inhouse) <- sDF_full_inhouse$Sentrix_ID

# # select betas of full inhouse and full cohort
sentrixs_full_inhouse <- intersect(sDF_full_inhouse$Sentrix_ID,rownames(betas))
sDF_full_inhouse <- sDF_full_inhouse[sDF_full_inhouse$Sentrix_ID %in% sentrixs_full_inhouse,]
Mset_full <- Mset[,sentrixs_full_inhouse]
betas_full <- betas[sentrixs_full_inhouse,]


message('saving clinical dataframe for full inhouse and full cohort in ', pDFclinical_full_inhouse)
saveRDS(sDF_full_inhouse, file = pDFclinical_full_inhouse)

message('saving clinical dataframe as cnv for for inhouse and full cohort in ', pDFclinical_full_inhouse_csv)
write.csv2(sDF_full_inhouse, file = pDFclinical_full_inhouse_csv)


message('saving Mset for full inhouse and full cohort in ', pMset_full)
saveRDS(Mset_full, file = pMset_full)

message('saving betas for full inhouse and full cohort in ', pbetas_full)
saveRDS(betas_full, file = pbetas_full)

# # select full cohort and 1p19q inhouse

isentrixs_full_gliomas <- which(sDF_full_inhouse$Cohort_ID %in% c('FullCohort', 'Selected1p19q')) 
sDF_full_gliomas <- sDF_full_inhouse[isentrixs_full_gliomas,]

message('saving clinical dataframe for full gliomas and full cohort in ', pDFclinical_full_gliomas)
saveRDS(sDF_full_gliomas, file = pDFclinical_full_gliomas)


# # select full cohort

isentrixs_cohort <- which(sDF_full_inhouse$Cohort_ID %in% c('FullCohort')) 
sDF_full_cohort <- sDF_full_inhouse[isentrixs_cohort,]

message('saving clinical dataframe for full cohort in ', pDFclinical_full_cohort)
saveRDS(sDF_full_cohort, file = pDFclinical_full_cohort)


message('saving clinical dataframe as cnv for full cohort in ', pDFclinical_full_cohort_csv)
write.csv2(sDF_full_cohort, file = pDFclinical_full_cohort_csv)


sink()
