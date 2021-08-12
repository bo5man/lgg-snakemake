
###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    # suppressMessages(library(minfi))
    suppressMessages(library(openxlsx))
    suppressMessages(library(dplyr))
    suppressMessages(library(ggplot2))
    suppressMessages(library(purrr))
    suppressMessages(library(RSpectra))
    suppressMessages(library(e1071))
    suppressMessages(library(randomForest))
    suppressMessages(library(randomForestSRC))
    suppressMessages(library(pROC))
} else{
    # library(minfi)
    library(openxlsx)
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(RSpectra)
    library(e1071)
    library(randomForest)
    library(randomForestSRC)
    library(pROC)
}

##################
# # # input file paths
pBetas_glass_treatment_related_620probes <-    snakemake@input[["betas_glass_treatment_related_620probes"]]
pBetas_glass_hypermodulator_342probes <-       snakemake@input[["betas_glass_hypermodulator_342probes"]]

pDFclinical_glass_treatment_related_620probes <-    snakemake@input[["DFclinical_glass_treatment_related_620probes"]]
pDFclinical_glass_hypermodulator_342probes <-       snakemake@input[["DFclinical_glass_hypermodulator_342probes"]]



# # # output file paths
# betas
pDFfeatures_glass_treatment_related_620probes <- snakemake@output[["DFfeatures_glass_treatment_related_620probes"]]
pforest_glass_treatment_related_620probes <-       snakemake@output[["forest_glass_treatment_related_620probes"]]
pforest_roc_glass_treatment_related_620probes <-   snakemake@output[["forest_roc_glass_treatment_related_620probes"]]

# # # parameters
dir_analysis_probes_rf <- snakemake@params[['dir_analysis_probes_rf']]
dir_analysis_probes <- snakemake@params[['dir_analysis_probes']]
dir_gCIMP <- snakemake@params[['dir_gCIMP']]
dir_glass <- snakemake@params[['dir_glass']]
dir_glass_treatment_related_probes <-   snakemake@params[['dir_glass_treatment_related_probes']]
dir_glass_hypermodulator_probes <-      snakemake@params[['dir_glass_hypermodulator_probes']]
#     dir_analysis_probes <- './../results/analysis_probes/'
#     dir_gCIMP <-                            paste0(dir_analysis_probes,'gCIMP_7probes/')
#     dir_glass <-                            paste0(dir_analysis_probes,'GLASS/')
#     dir_glass_treatment_related_probes <-   paste0(dir_glass,'treatment_related_620probes/')
#     dir_glass_hypermodulator_probes <-      paste0(dir_glass,'hypermodulator_342probes/')
#     pprobeset_gCIMP <- './../../LGG_Methylation/ProbeSet-gCIMP.xlsx'
#     pprobeset_GLASS <- './../../LGG_Methylation/ProbeSet-GLASS.xlsx'
#     probeset_GLASS_treatment_related_sheetname <-   'Treatment-related (N=620)'
#     probeset_GLASS_hypermodulator_sheetname <-      'IDHmut Hypermutator (N=342)'

if(!dir.exists(dir_analysis_probes_rf))(dir.create(dir_analysis_probes_rf))
if(!dir.exists(dir_gCIMP))(dir.create(dir_gCIMP))
if(!dir.exists(dir_glass))(dir.create(dir_glass))
if(!dir.exists(dir_glass_treatment_related_probes))(dir.create(dir_glass_treatment_related_probes))
if(!dir.exists(dir_glass_hypermodulator_probes))(dir.create(dir_glass_hypermodulator_probes))


log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################


#delete constant features
delete_constant_features <- function(myfeatures) {
  constant_features = c()
  for (j in 1:ncol(myfeatures)) {
    if (max(myfeatures[,j]) - min(myfeatures[,j]) == 0) {
      constant_features = c(constant_features,j)
    }
  }
  if (is.null(constant_features)) {
    return(myfeatures)
  } else {
    return(myfeatures[,-constant_features])
  }
}


#log transform of skewed features
deskew <- function(myfeatures, threshold) {
  myfeatures_deskewed = myfeatures
  for (j in 1:ncol(myfeatures)) {
    if (skewness(myfeatures[,j])>threshold & sum(myfeatures[,j]<=0)==0) { #skewed and positive
      for (i in 1:nrow(myfeatures)) {
        myfeatures_deskewed[i,j] = log(myfeatures[i,j])
      }
    }
  }
  return(myfeatures_deskewed)
}

# # # Read input
longshort <- readRDS(pDFclinical_glass_treatment_related_620probes)
betas <- readRDS(pBetas_glass_treatment_related_620probes)
betas <- betas[rownames(longshort),]

# # # select first 32000 betas with most spread 
# betas_most_spread <- betas[,order(-apply(betas,2,sd))[1:32000]]

# # final features for prediction!
features <- scale(deskew(delete_constant_features(betas), 1))
outcome <- factor(longshort$TypeSurvival, levels = c('long', 'short'))
DF_cohort <- data.frame(outcome, features)

message('saving dataframe containing outcome and features for full cohort in ', pDFfeatures_glass_treatment_related_620probes)
saveRDS(DF_cohort, file = pDFfeatures_glass_treatment_related_620probes)


# ##### RANDOM FOREST #####
# set.seed(0)
# 
# message('starting Random Forest at ', Sys.time())
# Forest <- rfsrc(outcome ~ .,data=DF_cohort,ntree=20000,var.used="all.trees",importance="TRUE",nodesize=2,seed=1)
# message('ending Random Forest at ', Sys.time())
# 
# message('saving Forest for full cohort in', pforest_full_cohort)
# saveRDS(Forest, file = pforest_full_cohort)
# 
# Forest_roc = pROC::roc(outcome~Forest$predicted.oob, levels=c('long','short'), direction="<")
# 
# message('saving ROC of Forest for full cohort in', pforest_roc_full_cohort)
# saveRDS(Forest_roc, file = pforest_roc_full_cohort)
# 
# plot(Forest$importance)
# 
# Forest_roc$auc #Out Of Bag performance
# Area under the curve: 0.6748
# Forest$importance #feature importances



sink()
sink()
