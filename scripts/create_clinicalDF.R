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
pbetas <- snakemake@input[["betas"]]
# pbetas <- './../results/betas/betas.RDS'

# # output file paths
pDFclinical_full_cohort <-      snakemake@output[["DFclinical_full_cohort"]]
pDFclinical_full_gliomas <-      snakemake@output[["DFclinical_full_gliomas"]]
pDFclinical_full_inhouse <-      snakemake@output[["DFclinical_full_inhouse"]]
# pDFclinical_full_gliomas <-      './../results/DFclinical_full_gliomas.RDS'

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
betas <- readRDS(pbetas)

# # # read clinical data files
inhouse_full <- read.xlsx(pinhouse, sep = "_")
cohort_full <- read.xlsx(pcohort, sep = '_')

# # # select inhouse data and full cohort data and create data frame sDF_inhouse
inhouse <- inhouse_full[,c('Sentrix_ID', 'Methylatie_array_uitkomst_schoon')]
# # # select cohort data
# cohort <- cohort_full[,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Shallow_Sequencing_Result', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]
cohort <- cohort_full[!is.na(cohort_full$SENTRIX_ID),c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Methylation_subclasses', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]

SelectedGliomas <- unique(str_subset(inhouse$Methylatie_array_uitkomst_schoon, 'Low Grade Glioma|1p/19q'))

iFullCohort <- seq_along(cohort$SENTRIX_ID)
iSelectedGliomas <- which(inhouse$Methylatie_array_uitkomst_schoon %in% SelectedGliomas)
iOtherSelectedGliomas <- setdiff(seq_along(inhouse$Methylatie_array_uitkomst_schoon), iSelectedGliomas)


sDF_full_inhouse <- data.frame(
                         ID = c(cohort$`Sample_ID_Long_vs_Short_Survivor_(LSS)`[iFullCohort],
                                 paste0('SG', seq_along(inhouse$Methylatie_array_uitkomst_schoon[iSelectedGliomas])),
                                 paste0('OG', seq_along(inhouse$Methylatie_array_uitkomst_schoon[iOtherSelectedGliomas]))
                                 ),
                         Cohort_ID = c(rep('FullCohort',length(cohort$SENTRIX_ID[iFullCohort])),
                                 rep('SelectedGliomas',length(inhouse$Sentrix_ID[iSelectedGliomas])),
                                 rep('OtherSelectedGliomas',length(inhouse$Sentrix_ID[iOtherSelectedGliomas]))
                                 ),
                         Sentrix_ID = c(cohort$SENTRIX_ID[iFullCohort],
                                 inhouse$Sentrix_ID[iSelectedGliomas],
                                 inhouse$Sentrix_ID[iOtherSelectedGliomas]
                                 ),
                         Type = c(cohort$Methylation_subclasses[iFullCohort],
                                 inhouse$Methylatie_array_uitkomst_schoon[iSelectedGliomas],
                                 inhouse$Methylatie_array_uitkomst_schoon[iOtherSelectedGliomas]
                                 ),
                         TypeSurvival = c(cohort$`long_or_short_(>100_months_-_<40_months)`[iFullCohort],
                                 inhouse$Methylatie_array_uitkomst_schoon[iSelectedGliomas],
                                 inhouse$Methylatie_array_uitkomst_schoon[iOtherSelectedGliomas]
                                 ),
                         SexType = c(cohort$`Sex_(M/F)`[iFullCohort],
                                 inhouse$Methylatie_array_uitkomst_schoon[iSelectedGliomas],
                                 inhouse$Methylatie_array_uitkomst_schoon[iOtherSelectedGliomas]
                                 )
)
rownames(sDF_full_inhouse) <- sDF_full_inhouse$Sentrix_ID

# # select betas of full inhouse and full cohort
sentrixs_full_inhouse <- intersect(sDF_full_inhouse$Sentrix_ID,rownames(betas))
sDF_full_inhouse <- sDF_full_inhouse[sDF_full_inhouse$Sentrix_ID %in% sentrixs_full_inhouse,]

message('saving clinical dataframe for full inhouse and full cohort in ', pDFclinical_full_inhouse)
saveRDS(sDF_full_inhouse, file = pDFclinical_full_inhouse)


# # select full cohort and gliomas inhouse

isentrixs_full_gliomas <- which(sDF_full_inhouse$Cohort_ID %in% c('FullCohort', 'SelectedGliomas')) 
sDF_full_gliomas <- sDF_full_inhouse[isentrixs_full_gliomas,]

message('saving clinical dataframe for full gliomas and full cohort in ', pDFclinical_full_gliomas)
saveRDS(sDF_full_gliomas, file = pDFclinical_full_gliomas)


# # select full cohort

isentrixs_cohort <- which(sDF_full_inhouse$Cohort_ID %in% c('FullCohort')) 
sDF_full_cohort <- sDF_full_inhouse[isentrixs_cohort,]

message('saving clinical dataframe for full cohort in ', pDFclinical_full_cohort)
saveRDS(sDF_full_cohort, file = pDFclinical_full_cohort)

sink()
