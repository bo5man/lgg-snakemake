###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(openxlsx))
    suppressMessages(library(dplyr))
    suppressMessages(library(gplots))
    suppressMessages(library(minfi))
    suppressMessages(library(purrr))
    suppressMessages(library(RSpectra))
    suppressMessages(library(Rtsne))
    suppressMessages(library(ggplot2))
} else{
    library(openxlsx)
    library(dplyr)
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
pinhouse <- snakemake@params[['inhouse']]
pcohort <- snakemake@params[['overview']]
# pinhouse <- './../../LGG_Methylation/Yongsoo VUMC methylation arrat samples sep 2020.xlsx'
# pcohort <- './../../LGG_Methylation/Overview_oligodendroglioma_samples_long_vs_short_PFS_2021-07-14.xlsx'

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

# read input
betas <- readRDS(pbetas)

# # # read clinical data files
inhouse_full <- read.xlsx(pinhouse, sep = "_")
cohort_full <- read.xlsx(pcohort, sep = '_')

# # # select cohort data and create data frame sDF
# cohort <- cohort_full[,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Shallow_Sequencing_Result', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]
cohort <- cohort_full[!is.na(cohort_full$SENTRIX_ID) ,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Methylation_subclasses', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]

sDF_full_cohort <- data.frame(
                         LSS = cohort$`Sample_ID_Long_vs_Short_Survivor_(LSS)`,
                         Sentrix_ID=cohort$SENTRIX_ID,
                         Type=cohort$Methylation_subclasses,
                         batch = rep('cohort',length(cohort$SENTRIX_ID)),
                         TypeSurvival=cohort$`long_or_short_(>100_months_-_<40_months)`,
                         SexType=cohort$`Sex_(M/F)`
)
rownames(sDF_full_cohort) <- sDF_full_cohort$Sentrix_ID

i <- intersect(rownames(sDF_full_cohort),rownames(betas))
sDF_full_cohort <- sDF_full_cohort[i,]

message('saving clinical dataframe for full cohort in ', pDFclinical_full_cohort)
saveRDS(sDF_full_cohort, file = pDFclinical_full_cohort)



# # # select inhouse gliomas and create data frame sDF_gliomas
inhouse <- inhouse_full[,c('Sentrix_ID', 'Methylatie_array_uitkomst_schoon')]
# # # select cohort data
# cohort <- cohort_full[,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Shallow_Sequencing_Result', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]
cohort <- cohort_full[!is.na(cohort_full$SENTRIX_ID),c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Methylation_subclasses', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]


# select relevant Gliomas
TypesGliomas <- unique(grep('Low Grade Glioma|1p/19q', inhouse$Methylatie_array_uitkomst_schoon, value= T))
inhouse_gliomas <- filter(inhouse, Methylatie_array_uitkomst_schoon %in% TypesGliomas)

sDF_full_gliomas <- data.frame(
                         LSS = c(cohort$`Sample_ID_Long_vs_Short_Survivor_(LSS)`,inhouse_gliomas$Methylatie_array_uitkomst_schoon),
                         Sentrix_ID=c(cohort$SENTRIX_ID,inhouse_gliomas$Sentrix_ID),
                         Type=c(cohort$Methylation_subclasses,inhouse_gliomas$Methylatie_array_uitkomst_schoon),
                         TypeSurvival=c(cohort$`long_or_short_(>100_months_-_<40_months)`,inhouse_gliomas$Methylatie_array_uitkomst_schoon),
                         SexType=c(cohort$`Sex_(M/F)`,inhouse_gliomas$Methylatie_array_uitkomst_schoon)
)
rownames(sDF_full_gliomas) <- sDF_full_gliomas$Sentrix_ID

# # select betas of gliomas
i <- intersect(rownames(sDF_full_gliomas),rownames(betas))
sDF_full_gliomas <- sDF_full_gliomas[i,]

message('saving clinical dataframe for Gliomas and full cohort in ', pDFclinical_full_gliomas)
saveRDS(sDF_full_gliomas, file = pDFclinical_full_gliomas)



# # # select inhouse data and full cohort data and create data frame sDF_inhouse
inhouse <- inhouse_full[,c('Sentrix_ID', 'Methylatie_array_uitkomst_schoon')]
# # # select cohort data
# cohort <- cohort_full[,c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Shallow_Sequencing_Result', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]
cohort <- cohort_full[!is.na(cohort_full$SENTRIX_ID),c('Sample_ID_Long_vs_Short_Survivor_(LSS)','SENTRIX_ID','Sufficient_DNA_and_paid_for', 'Methylation_subclasses', 'long_or_short_(>100_months_-_<40_months)', 'Sex_(M/F)')]

sDF_full_inhouse <- data.frame(
                         LSS = c(cohort$`Sample_ID_Long_vs_Short_Survivor_(LSS)`,inhouse$Methylatie_array_uitkomst_schoon),
                         Sentrix_ID=c(cohort$SENTRIX_ID,inhouse$Sentrix_ID),
                         Type=c(cohort$Methylation_subclasses,inhouse$Methylatie_array_uitkomst_schoon),
                         TypeSurvival=c(cohort$`long_or_short_(>100_months_-_<40_months)`,inhouse$Methylatie_array_uitkomst_schoon),
                         SexType=c(cohort$`Sex_(M/F)`,inhouse$Methylatie_array_uitkomst_schoon)
)
rownames(sDF_full_inhouse) <- sDF_full_inhouse$Sentrix_ID

# # select betas of inhouse and full cohort
i <- intersect(rownames(sDF_full_inhouse),rownames(betas))
sDF_full_inhouse <- sDF_full_inhouse[i,]

message('saving clinical dataframe for inhouse and full cohort in ', pDFclinical_full_inhouse)
saveRDS(sDF_full_inhouse, file = pDFclinical_full_inhouse)



sink()
sink()
