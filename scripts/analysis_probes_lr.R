
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
}

##################
# # # input file paths
pBetas_gCIMP <- snakemake@input[["betas_gCIMP"]]
pBetas_gCIMP <- snakemake@input[["betas_GLASS_treatment_related"]]
#pBetas_gCIMP_thresholded <- snakemake@input[["betas_gCIMP_thresholded"]]
# pBetas_gCIMP <- './../results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS'
# pBetas_gCIMP_thresholded <- './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_thresholded.RDS'

pDFclinical_cohort <-        snakemake@input[["DFclinical_cohort"]]
# pDFclinical_cohort <- './../results/DFclinical_cohort.RDS'

# # # output file paths
# betas
pglm_stats_gCIMP <- snakemake@output[["glm_stats_gCIMP"]]
# pglm_stats_gCIMP_thresholded <- snakemake@output[["glm_stats_gCIMP_thresholded"]]
pDFfeatures_gCIMP <- snakemake@output[["DFfeatures_gCIMP"]]
# pDFfeatures_gCIMP_thresholded <- snakemake@output[["DFfeatures_gCIMP_thresholded"]]
# pglm_stats_gCIMP <- "./../results/analysis_probes_lr/glm_stats_gCIMP.RDS"

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




# # # Read input
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


# # # final thresholded features for prediction!
# features_thresholded <- scale(deskew(delete_constant_features(betas_thresholded), 1))
# outcome <- factor(longshort$TypeSurvival, levels = c('long', 'short'))
# DF_cohort_thresholded <- data.frame(outcome, features_thresholded)
# 
# message('saving dataframe containing outcome and thresholded features for full cohort with gCIMP probes for linear regression in  ', pDFfeatures_gCIMP_thresholded)
# saveRDS(DF_cohort_thresholded, file = pDFfeatures_gCIMP_thresholded)


# # # set seed for reproducability
set.seed(0)


# # # General linear model for gCIMP features
myalpha = 0 #alpha=0 ridge, alpha=1 LASSO
penaltyfactor = rep(1,ncol(features)) #1 is penalized, 0 is unpenalized
train_percentage = 0.8 
nrepeats = 250

aucs = numeric(nrepeats)
myrocs = list()
model_trains = list()
prob_tests = list()
indices0 = (1:length(outcome))[outcome=='short']
indices1 = (1:length(outcome))[outcome=='long']
for (k in 1:nrepeats) {
  ind0 = sample(indices0)
  ind1 = sample(indices1)
  indices_train = c(ind0[1:round(length(ind0)*train_percentage)],
                    ind1[1:round(length(ind1)*train_percentage)])
  indices_test = c(ind0[(round(length(ind0)*train_percentage)+1):length(ind0)],
                   ind1[(round(length(ind1)*train_percentage)+1):length(ind1)])
  net = cv.glmnet(features[indices_train,], outcome[indices_train], alpha=myalpha, penalty.factor = penaltyfactor, nfolds=5, type.measure='class', family="binomial")
  # plot(net)
  best_lambda = net$lambda.min
  model_trains[[k]] = glmnet(features[indices_train,], outcome[indices_train], alpha=myalpha, penalty.factor = penaltyfactor, family="binomial", lambda=best_lambda)
  prob_tests[[k]] <- predict(model_trains[[k]], newx=features[indices_test,], type="response", type.measure='class')
  myrocs[[k]] = pROC::roc(outcome[indices_test]~as.vector(prob_tests[[k]]),
                    #                     levels=c('short','long'),
                    direction="<", quiet=TRUE)
  aucs[k] = myrocs[[k]]$auc
}

# message('The gCIMP features model: the mean area under the curve over ',nrepeats,' repeats is ', mean(aucs))
# message('saving the trained generalized linear model, test probabilities and roc for gCIMP features in ', pglm_stats_gCIMP)
# glm_stats_gCIMP = list('model_trains' = model_trains, 'prob_tests' = prob_tests, 'myrocs' = myrocs)
# saveRDS(glm_stats_gCIMP, file = pglm_stats_gCIMP)


# # # sample indices (train_percentage*nrow(outcome))x nrepeats)-matrix


cvnets = list()
model_trains = list()
prediction_validations = list()


indices_train_times <- createDataPartition(outcome,times= nrepeats, p=train_percentage, list=F)
for (k in 1:nrepeats){
    indices_train <- indices_train_times[,k]
    training_expression <- features[indices_train,]
    training_phenotype <- outcome[indices_train]
    validation_expression <- features[-indices_train,]
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
    # as.table(cm)
    # write.xlsx(as.table(cm),file=pglm_stats_gCIMP, sheetName='confusionMatrix')
    # # statistics like Sensitivity, Specificity etc
    # as.matrix(cm$byClass, what = "classes")
    # write.xlsx(as.matrix(cm$byClass, what = "classes"),file=pglm_stats_gCIMP, sheetName='statistics')
    
}

# message('saving the trained generalized linear model, test probabilities and roc for gCIMP features in ', pglm_stats_gCIMP)
# glm_stats_gCIMP = list('model_trains' = model_trains, 'prob_tests' = prob_tests, 'myrocs' = myrocs)
# saveRDS(glm_stats_gCIMP, file = pglm_stats_gCIMP)




sink()
