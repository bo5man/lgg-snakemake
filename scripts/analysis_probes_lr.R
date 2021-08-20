
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
    suppressMessages(library(glmnet))
    suppressMessages(library(randomForest))
    suppressMessages(library(randomForestSRC))
    suppressMessages(library(pROC))
    suppressMessages(library(caret))
} else{
    # library(minfi)
    library(openxlsx)
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(RSpectra)
    library(e1071)
    library(glmnet)
    library(randomForest)
    library(randomForestSRC)
    library(pROC)
    library(caret)
}

##################
# # # input file paths
pBetas_gCIMP <- snakemake@input[["betas_gCIMP"]]
pBetas_glass_treatment_related <- snakemake@input[["betas_glass_treatment_related_620probes"]]
#pBetas_gCIMP <- snakemake@input[["betas_GLASS_treatment_related"]]
#pBetas_gCIMP_thresholded <- snakemake@input[["betas_gCIMP_thresholded"]]
# pBetas_gCIMP <- './../../output-lgg-test/results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS'
# pBetas_gCIMP_thresholded <- './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_thresholded.RDS'

pDFclinical_cohort <-        snakemake@input[["DFclinical_cohort"]]
# pDFclinical_cohort <- './../results/DFclinical_cohort.RDS'

# # # output file paths
# betas
pglm_fullstats_gCIMP <- snakemake@output[["glm_fullstats_gCIMP"]]
pglm_stats_gCIMP <- snakemake@output[["glm_stats_gCIMP"]]
pDFfeatures_gCIMP <- snakemake@output[["DFfeatures_gCIMP"]]
pglm_fullstats_glass_treatment_related <- snakemake@output[["glm_fullstats_glass_treatment_related"]]
pglm_stats_glass_treatment_related <- snakemake@output[["glm_stats_glass_treatment_related"]]
pDFfeatures_glass_treatment_related <- snakemake@output[["DFfeatures_glass_treatment_related"]]

# # # parameters
dir_analysis_probes_lr <- snakemake@params[['dir_analysis_probes_lr']]
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

if(!dir.exists(dir_analysis_probes_lr))(dir.create(dir_analysis_probes_lr))
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




# # # Read input gCIMP data
longshort <- readRDS(pDFclinical_cohort)
betas <- readRDS(pBetas_gCIMP)
betas <- betas[rownames(longshort),]
# longshort = readRDS('./../results/DFclinical_cohort.RDS')
# betas = readRDS('./../results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS')
# betas_thresholded <- readRDS(pBetas_gCIMP_thresholded)
# betas_thresholded <- betas_thresholded[rownames(longshort),]

warning('removing sample 205061430022_R06C01(LSS15) from cohort due to high grade astrocytoma, skewing the features, see heatmaps')
longshort <- longshort[setdiff(rownames(longshort),'205061430022_R06C01'),]
betas <- betas[setdiff(rownames(betas),'205061430022_R06C01'),]

# # final features for prediction!
features <- scale(deskew(delete_constant_features(betas), 1))
outcome <- factor(longshort$TypeSurvival, levels = c('long', 'short'))
DF_cohort <- data.frame(outcome, features)

message('saving dataframe containing outcome and features for cohort with gCIMP probes for linear regression in ', pDFfeatures_gCIMP)
saveRDS(DF_cohort, file = pDFfeatures_gCIMP)


# # # set seed for reproducability
set.seed(0)


# # # General linear model for gCIMP features
myalpha = 0 #alpha=0 ridge, alpha=1 LASSO
penaltyfactor = rep(1,ncol(features)) #1 is penalized, 0 is unpenalized
train_percentage = 0.8 
nrepeats = 250

# # # sample indices (train_percentage*nrow(outcome))x nrepeats)-matrix
cvnets = list()
model_trains = list()
prediction_validations = list()
cms <- list()
cms_stats <- matrix(,nrow=nrepeats,ncol=18)
colnames(cms_stats) <- c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyNull', 'AccuracyPValue', 'McnemarPValue', 'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value', 'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate', 'Detection Prevalence', 'Balanced Accurary')
rownames(cms_stats) <- paste0('repeat',1:nrepeats)

indices_train_times <- createDataPartition(outcome,times= nrepeats, p=train_percentage, list=F)
for (k in 1:nrepeats){
    indices_train <- indices_train_times[,k]
    training_expression <- as.matrix(features[indices_train,])
    training_phenotype <- outcome[indices_train]
    validation_expression <- as.matrix(features[-indices_train,])
    validation_phenotype <- outcome[-indices_train]
    # # cross validation
    cvnet <- cv.glmnet(training_expression, training_phenotype, alpha=myalpha, penalty.factor = penaltyfactor, nfolds=5, type.measure='class', family="binomial")
    cvnets[[k]] <- cvnet
    # plot(cvnet)
    best_lambda <- cvnet$lambda.min
    # # actual model training
    model_train <- glmnet(training_expression, training_phenotype, alpha=myalpha, penalty.factor = penaltyfactor, family="binomial", lambda=best_lambda)
    model_trains[[k]] <- model_train
    # # prediction of classes for validation set
    prediction_validation <- predict(model_train, newx=validation_expression, type="class", type.measure='class')
    prediction_validations[[k]] <- prediction_validation
    prediction_validation <- factor(as.vector(prediction_validation),levels=levels(outcome))
    # # confusion Matrix
    cm <- confusionMatrix(prediction_validation,validation_phenotype)
    cms[[k]] <- cm
    cms_stats[k,] <- c(cm$overall,cm$byClass)
}

message('saving all details of confusion matrix in ', pglm_fullstats_gCIMP)
saveRDS(cms, file = pglm_fullstats_gCIMP)
cms_stats_average <- as.matrix(colMeans(cms_stats,na.rm = T))
colnames(cms_stats_average) <- 'statistic'

message('saving statistics and average statistics of confusion matrix in ', pglm_stats_gCIMP)
wb <- createWorkbook()
addWorksheet(wb, "statistics")
addWorksheet(wb, "average_statistics")

writeData(wb, "statistics", cms_stats, startRow = 1, startCol = 1, rowNames = T)
writeData(wb, "average_statistics", cms_stats_average, startRow = 1, startCol = 1, rowNames = T)

saveWorkbook(wb, file = pglm_stats_gCIMP, overwrite = TRUE)




# # # Read input glass treatment related data
longshort <- readRDS(pDFclinical_cohort)
betas <- readRDS(pBetas_glass_treatment_related)
betas <- betas[rownames(longshort),]
# longshort = readRDS('./../results/DFclinical_cohort.RDS')
# betas = readRDS('./../results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS')
# betas_thresholded <- readRDS(pBetas_gCIMP_thresholded)
# betas_thresholded <- betas_thresholded[rownames(longshort),]

warning('removing sample 205061430022_R06C01(LSS15) from cohort due to high grade astrocytoma, skewing the features, see heatmaps')
longshort <- longshort[setdiff(rownames(longshort),'205061430022_R06C01'),]
betas <- betas[setdiff(rownames(betas),'205061430022_R06C01'),]

# # final features for prediction!
features <- scale(deskew(delete_constant_features(betas), 1))
outcome <- factor(longshort$TypeSurvival, levels = c('long', 'short'))
DF_cohort <- data.frame(outcome, features)

message('saving dataframe containing outcome and features for cohort with glass treatment related probes for linear regression in ', pDFfeatures_glass_treatment_related)
saveRDS(DF_cohort, file = pDFfeatures_glass_treatment_related)


# # # set seed for reproducability
set.seed(0)


# # # General linear model for gCIMP features
myalpha = 0 #alpha=0 ridge, alpha=1 LASSO
penaltyfactor = rep(1,ncol(features)) #1 is penalized, 0 is unpenalized
train_percentage = 0.8 
nrepeats = 250

# # # sample indices (train_percentage*nrow(outcome))x nrepeats)-matrix
cvnets = list()
model_trains = list()
prediction_validations = list()
cms <- list()
cms_stats <- matrix(,nrow=nrepeats,ncol=18)
colnames(cms_stats) <- c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyNull', 'AccuracyPValue', 'McnemarPValue', 'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value', 'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate', 'Detection Prevalence', 'Balanced Accurary')
rownames(cms_stats) <- paste0('repeat',1:nrepeats)

indices_train_times <- createDataPartition(outcome,times= nrepeats, p=train_percentage, list=F)
for (k in 1:nrepeats){
    indices_train <- indices_train_times[,k]
    training_expression <- as.matrix(features[indices_train,])
    training_phenotype <- outcome[indices_train]
    validation_expression <- as.matrix(features[-indices_train,])
    validation_phenotype <- outcome[-indices_train]
    # # cross validation
    cvnet <- cv.glmnet(training_expression, training_phenotype, alpha=myalpha, penalty.factor = penaltyfactor, nfolds=5, type.measure='class', family="binomial")
    cvnets[[k]] <- cvnet
    # plot(cvnet)
    best_lambda <- cvnet$lambda.min
    # # actual model training
    model_train <- glmnet(training_expression, training_phenotype, alpha=myalpha, penalty.factor = penaltyfactor, family="binomial", lambda=best_lambda)
    model_trains[[k]] <- model_train
    # # prediction of classes for validation set
    prediction_validation <- predict(model_train, newx=validation_expression, type="class", type.measure='class')
    prediction_validations[[k]] <- prediction_validation
    prediction_validation <- factor(as.vector(prediction_validation),levels=levels(outcome))
    # # confusion Matrix
    cm <- confusionMatrix(prediction_validation,validation_phenotype)
    cms[[k]] <- cm
    cms_stats[k,] <- c(cm$overall,cm$byClass)
}

message('saving all details of confusion matrix in ', pglm_fullstats_glass_treatment_related)
saveRDS(cms, file = pglm_fullstats_glass_treatment_related)
cms_stats_average <- as.matrix(colMeans(cms_stats,na.rm = T))
colnames(cms_stats_average) <- 'statistic'

message('saving statistics and average statistics of confusion matrix in ', pglm_stats_glass_treatment_related)
wb <- createWorkbook()
addWorksheet(wb, "statistics")
addWorksheet(wb, "average_statistics")

writeData(wb, "statistics", cms_stats, startRow = 1, startCol = 1, rowNames = T)
writeData(wb, "average_statistics", cms_stats_average, startRow = 1, startCol = 1, rowNames = T)

saveWorkbook(wb, file = pglm_stats_glass_treatment_related, overwrite = TRUE)


sink()
