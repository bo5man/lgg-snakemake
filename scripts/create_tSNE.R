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
    suppressMessages(library(ggfortify))
} else{
    library(openxlsx)
    library(dplyr)
    library(gplots)
    library(minfi)
    library(purrr)
    library(RSpectra)
    library(Rtsne)
    library(ggplot2)
    library(ggfortify)
}

##################
# # input file paths
pbetas_full <- snakemake@input[["betas_full"]]
# pBetas_cohort<- snakemake@input[["betas_cohort"]]
# pBetas_gliomas <- snakemake@input[["betas_gliomas"]]
pBetas_inhouse <- snakemake@input[["betas_inhouse"]]
# pDFclinical_cohort <-  snakemake@input[["DFclinical_cohort"]]
# pDFclinical_gliomas <-      snakemake@input[["DFclinical_gliomas"]]
pDFclinical_inhouse <-      snakemake@input[["DFclinical_inhouse"]]
pDFclinical_full_inhouse <-      snakemake@input[["DFclinical_full_inhouse"]]
# qcMset <- snakemake@input[['qcMset']]
# pbetas <- './../results/betas/betas.RDS'
# pqcplotMset_cohort <- './../results/Mset/qcplotMset_cohort.png'

# # output file paths
# ppca_cohort <- snakemake@output[["pca_cohort"]]
# ppca_gliomas <- snakemake@output[["pca_gliomas"]]
ppca_inhouse <- snakemake@output[["pca_inhouse"]]
# ptsne_cohort <- snakemake@output[["tsne_cohort"]]
# ptsne_gliomas <- snakemake@output[["tsne_gliomas"]]
ptsne_inhouse <- snakemake@output[["tsne_inhouse"]]
ptsne_full_inhouse <- snakemake@output[["tsne_full_inhouse"]]
# ptsneplot_cohort_Type <- snakemake@output[["tsneplot_cohort_Type"]]
# ptsneplot_cohort_TypeSurvival <- snakemake@output[["tsneplot_cohort_TypeSurvival"]]
# ptsneplot_cohort_ID <- snakemake@output[["tsneplot_cohort_ID"]]
# ptsneplot_gliomas_Type <- snakemake@output[["tsneplot_gliomas_Type"]]
# ptsneplot_gliomas_TypeSurvival <- snakemake@output[["tsneplot_gliomas_TypeSurvival"]]
# ptsneplot_gliomas_ID <- snakemake@output[["tsneplot_gliomas_ID"]]
ptsneplot_inhouse_Type <-           snakemake@output[["tsneplot_inhouse_Type"]]
ptsneplot_inhouse_TypeSurvival <-   snakemake@output[["tsneplot_inhouse_TypeSurvival"]]
ptsneplot_inhouse_ID <-            snakemake@output[["tsneplot_inhouse_ID"]]
ptsneplot_inhouse_BadSample <- snakemake@output[["tsneplot_inhouse_BadSample"]]


# # parameters

dir_tsne <- snakemake@params[['dir_tsne']]
# dir_perplexity <- snakemake@params[['dir_perplexity']]
dir_rds <- snakemake@params[['dir_rds']]
dir_pca <- snakemake@params[['dir_pca']]
# dir_mset <- snakemake@params[['dir_mset']]
# dir_pca <- './../results/pca/'
if(!dir.exists(dir_pca))(dir.create(dir_pca))
if(!dir.exists(dir_tsne))(dir.create(dir_tsne))
if(!dir.exists(dir_rds))(dir.create(dir_rds))
# if(!dir.exists(dir_mset))(dir.create(dir_mset))

# pinhouse <- snakemake@params[['inhouse']]
# pcohort <- snakemake@params[['overview']]
# # pinhouse <- './../../LGG_Methylation/Yongsoo VUMC methylation arrat samples sep 2020.xlsx'
# # pcohort <- './../../LGG_Methylation/Overview_oligodendroglioma_samples_long_vs_short_PFS_2021-06-10.xlsx'

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# # # read input
# betas_cohort <- readRDS(pBetas_cohort)
# betas_gliomas <- readRDS(pBetas_gliomas)
betas_inhouse <- readRDS(pBetas_inhouse)

# sDF_cohort <- readRDS(pDFclinical_cohort)
# sDF_gliomas <- readRDS(pDFclinical_gliomas)
sDF_inhouse <- readRDS(pDFclinical_inhouse)


# # sanity check with bad samples
betas_full_inhouse <- readRDS(pbetas_full)
sDF_full_inhouse <- readRDS(pDFclinical_full_inhouse)

all(rownames(betas_full_inhouse)==rownames(sDF_full_inhouse))
# [1] TRUE
all(rownames(betas_inhouse)==rownames(sDF_inhouse))
# [1] TRUE

# # obtain sentrix and indici of bad samples to filter out
SentrixBad <- setdiff(rownames(betas_full_inhouse),rownames(betas_inhouse))
iSentrixBad <- which(rownames(betas_full_inhouse) %in% SentrixBad)

message('temporarily overwriting betas_inhouse with betas_full_inhouse and sDF_inhouse with sDF_full_inhouse')
betas_inhouse <- betas_full_inhouse
sDF_inhouse <- sDF_full_inhouse

# # Set seed for reproducability
set.seed(0)


# # calculate pcas for inhouse and cohort
kpca <- min(nrow(betas_inhouse),ncol(betas_inhouse))-1
pca_inhouse <- prcomp_svds(betas_inhouse,k=kpca)

message('saving pca for inhouse and cohort in ', ppca_inhouse)
saveRDS(pca_inhouse, file = ppca_inhouse)


# # calculate tSNE for inhouse and cohort

# default perplexity is 30, and should satisfy 3*perp < nrow(X)-1
perp = min(30,(nrow(betas_inhouse)-4)/3)

tsne_inhouse <- Rtsne(pca_inhouse$x,pca=F,max_iter=5500,theta=0,verbose=T,perplexity = perp)
message('saving tsne of inhouse and cohort in ', ptsne_inhouse)
saveRDS(tsne_inhouse, file = ptsne_inhouse)

# library(ggplot2)
# sDF_inhouse = readRDS('./../../output-lgg-test/results/DFclinical_full_inhouse.RDS')
# tsne_inhouse=readRDS('./../../output-lgg-test/results/tSNE/RDS/tsne_inhouse.RDS')
# perp=1;kpca=1;

sDF_inhouse$X1 <- tsne_inhouse$Y[,1]
sDF_inhouse$X2 <- tsne_inhouse$Y[,2]


# plot: FullCohort vs Selected1p19q vs OtherGliomas
# plot: FullCohort vs Selected1p19q vs OtherGliomas, which samples are bad in full cohort. labeled SentrixBad red, good green pan to cnvs of bad and good
# plot: FullCohort vs Selected1p19q vs OtherGliomas, where do cohort 1p19q lay compared to inhouse 1p19q.

# plot: FullCohort vs Selected1p19q vs OtherGliomas Cohort_ID
message('saving ptsneplot_inhouse_ID in ', ptsneplot_inhouse_ID)
png(file = ptsneplot_inhouse_ID, width = 1440, height = 1280) #width = 1846,height = 991
ggplot(sDF_inhouse, aes(x=X1, y=X2, shape=Cohort_ID, colour=Cohort_ID)) +
    theme_bw() + 
    geom_point() +
    geom_text(aes(x=X1,y=X2,label=ID), show.legend = F, nudge_y=0.25, check_overlap=T) +
    guides(shape=guide_legend(ncol=3)) +
    guides(label="none") +
    theme(legend.position = "bottom") +
    ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE)-plot of the principal component analysis(pca) of all methylation probes of the full cohort(LSS), selected 1p/19q inhouse sample(1p19qI) and other related inhouse glioma samples(OG)."),
            subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
    )
dev.off()

# plot: FullCohort vs Selected1p19q vs OtherGliomas Type
message('saving tsneplot_inhouse_Type in ', ptsneplot_inhouse_Type)
png(file = ptsneplot_inhouse_Type, width = 1440, height = 1280) #width = 1846,height = 991
ggplot(sDF_inhouse, aes(x=X1, y=X2)) +
    theme_bw() + 
    geom_point(aes(x=X1,y=X2, shape=Cohort_ID), show.legend=F) +
    geom_text(aes(x=X1,y=X2, label=Type, col=Type), show.legend = T, nudge_y=0.25, check_overlap=T) +
    guides(col=guide_legend(ncol=4)) +
    guides(shape="none") +
    #     guides(shape=guide_legend(ncol=1)) +
    theme(legend.position = "bottom") +
    ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE)-plot of the principal component analysis(pca) of all methylation probes of the full cohort(LSS), selected 1p/19q inhouse sample(1p19qI) and other related inhouse glioma samples(OG)."),
            subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
    )
dev.off()

# plot: FullCohort vs Selected1p19q vs OtherGliomas, which samples are bad in full cohort. labeled SentrixBad red, good green pan to cnvs of bad and good
message('saving tsneplot_inhouse_BadSample in ', ptsneplot_inhouse_BadSample)
png(file = ptsneplot_inhouse_BadSample, width = 1440, height = 1280) #width = 1846,height = 991
ggplot(sDF_inhouse, aes(x=X1, y=X2, shape=Cohort_ID, colour=Cohort_ID)) +
    theme_bw() + 
    geom_point()  +
    geom_text(data = sDF_inhouse[sDF_inhouse$Sentrix_ID %in% SentrixBad,],aes(x=X1,y=X2,label=Type),colour="black",nudge_y=0.25) +
    guides(shape=guide_legend(ncol=3)) +
    guides(label="none") +
    theme(legend.position = "bottom") +
    ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE)-plot of the principal component analysis(pca) of all methylation probes of the full cohort, selected 1p/19q inhouse sample and other related inhouse glioma samples.\nLabeled samples depict samples from the full cohort that did not receive the 1p/19q subclass by the methylation classifier."),
            subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
    )
dev.off()


# plot: FullCohort vs Selected1p19q, where do cohort 1p19q lay compared to inhouse 1p19q.
message('saving tsneplot_inhouse_TypeSurvival in ', ptsneplot_inhouse_TypeSurvival)
png(file = ptsneplot_inhouse_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
ggplot(sDF_inhouse[!sDF_inhouse$Sentrix_ID %in% SentrixBad & sDF_inhouse$Cohort_ID != 'OtherSelectedGliomas',],aes(x=X1, y=X2, shape=Cohort_ID))  +
    theme_bw() + 
    geom_point()  +
    geom_text(data = sDF_inhouse[!sDF_inhouse$Sentrix_ID %in% SentrixBad & sDF_inhouse$Cohort_ID ==  'FullCohort',], aes(x=X1,y=X2, label=TypeSurvival, col=TypeSurvival),nudge_y=0.25, check_overlap=T) +
    guides(shape=guide_legend(ncol=2)) +
    guides(col=guide_legend(ncol=2)) +
    theme(legend.position = "bottom") +
    ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE)-plot of the principal component analysis(pca) of all methylation probes of the full cohort, selected 1p/19q inhouse sample and other related inhouse glioma samples.\nShown are the selected 1p/19q inhouse samples and the 1p/19q samples of the cohort, labeled by their survival."),
            subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
    )
dev.off()



# # # # # # # #

# DEPRECATED

# # # # # # # #


# message('saving tsneplot_inhouse_Type in ', ptsneplot_inhouse_Type)
# png(file = ptsneplot_inhouse_Type, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_inhouse, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     # geom_text(data = sDF_inhouse,aes(x=X1,y=X2,label=Type)) +
#     geom_text(data = sDF_inhouse,aes(x=X1,y=X2,label=ID, col=Type), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=2)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and related inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# 
# # # Example of showing only the Cohort
# 
# 
# 
# 
# message('saving tsneplot_inhouse_TypeSurvival in ', ptsneplot_inhouse_TypeSurvival)
# png(file = ptsneplot_inhouse_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_inhouse, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     geom_text(data = sDF_inhouse,aes(x=X1,y=X2,label=TypeSurvival), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=2)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and related inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# message('saving tsneplot_inhouse_ID in ', ptsneplot_inhouse_ID)
# png(file = ptsneplot_inhouse_ID, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_inhouse, aes(x=X1, y=X2, col=TypeSurvival)) +
#     geom_point()  +
#     geom_text(data = sDF_inhouse,aes(x=X1,y=X2,label=ID), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and related inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# 
# # # calculate pcas for cohort
# kpca <- min(nrow(betas_cohort),ncol(betas_cohort))-1
# # all(rownames(betas_cohort) == rownames(sDF_cohort))
# rownames(sDF_cohort) <- sDF_cohort$ID 
# rownames(betas_cohort) <- sDF_cohort$ID 
# pca_cohort <- prcomp_svds(betas_cohort,k=kpca)
# 
# message('saving pca for cohort in ', ppca_cohort)
# saveRDS(pca_cohort, file = ppca_cohort)
# 
# 
# # # # show principle components
# # 
# # message('saving pca of cohort in ', ppcaplot_cohort_TypeSurvival)
# # png(file = ppcaplot_cohort_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
# # p1 <- autoplot(pca_cohort, label = TRUE, data = sDF_cohort , colour = 'TypeSurvival') +
# #     labs(title = 'Principal Component Analysis of 850k probes for 19 cohort samples of 1p/19q oligodendroglioma')
# # dev.off()
# 
# # # calculate tSNE for cohort
# 
# # default perplexity is 30, and should satisfy 3*perp < nrow(X)-1
# perp = min(30,(nrow(betas_cohort)-4)/3)
# 
# tsne_cohort <- Rtsne(pca_cohort$x,pca=F,max_iter=5500,theta=0,verbose=T,perplexity = perp)
# message('saving tsne of cohort in ', ptsne_cohort)
# saveRDS(tsne_cohort, file = ptsne_cohort)
# 
# sDF_cohort$X1 <- tsne_cohort$Y[,1]
# sDF_cohort$X2 <- tsne_cohort$Y[,2]
# 
# message('saving tsneplot_cohort_Type in ', ptsneplot_cohort_Type)
# png(file = ptsneplot_cohort_Type, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_cohort, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     geom_text(data = sDF_cohort,aes(x=X1,y=X2,label=ID, col=Type), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# message('saving tsneplot_cohort_TypeSurvival in ', ptsneplot_cohort_TypeSurvival)
# png(file = ptsneplot_cohort_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_cohort, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     geom_text(data = sDF_cohort,aes(x=X1,y=X2,label=TypeSurvival), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# message('saving tsneplot_cohort_ID in ', ptsneplot_cohort_ID)
# png(file = ptsneplot_cohort_ID, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_cohort, aes(x=X1, y=X2, col=TypeSurvival)) +
#     geom_point()  +
#     geom_text(data = sDF_cohort,aes(x=X1,y=X2,label=ID), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# 
# 
# # # calculate pcas for gliomas
# kpca <- min(nrow(betas_gliomas),ncol(betas_gliomas))-1
# 
# # all(rownames(betas_cohort) == rownames(sDF_cohort))
# # rownames(sDF_gliomas) <- sDF_gliomas$ID
# pca_gliomas <- prcomp_svds(betas_gliomas,k=kpca)
# 
# message('saving pca for gliomas in ', ppca_gliomas)
# saveRDS(pca_gliomas, file = ppca_gliomas)
# 
# # # # show principle components
# # 
# # message('saving pca of gliomas in ', ppcaplot_gliomas_TypeSurvival)
# # png(file = ppcaplot_gliomas_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
# # p2 <- autoplot(pca_gliomas, label = TRUE, data = sDF_gliomas , colour = 'TypeSurvival') +
# #     labs(title = 'Principal Component Analysis of 850k probes for 19 cohort samples of 1p/19q oligodendroglioma/n and related glioma inhouse samples')
# # dev.off()
# 
# # # calculate tSNE for gliomas
# 
# # default perplexity is 30, and should satisfy 3*perp < nrow(X)-1
# perp = min(30,(nrow(betas_gliomas)-4)/3)
# 
# tsne_gliomas <- Rtsne(pca_gliomas$x,pca=F,max_iter=5500,theta=0,verbose=T,perplexity = perp)
# message('saving tsne of gliomas and cohort in ', ptsne_gliomas)
# saveRDS(tsne_gliomas, file = ptsne_gliomas)
# 
# sDF_gliomas$X1 <- tsne_gliomas$Y[,1]
# sDF_gliomas$X2 <- tsne_gliomas$Y[,2]
# 
# message('saving tsneplot_gliomas_Type in ', ptsneplot_gliomas_Type)
# png(file = ptsneplot_gliomas_Type, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_gliomas, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     geom_text(data = sDF_gliomas,aes(x=X1,y=X2,label=ID, col=Type), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and selected 1p/19q inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# message('saving tsneplot_gliomas_TypeSurvival in ', ptsneplot_gliomas_TypeSurvival)
# png(file = ptsneplot_gliomas_TypeSurvival, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_gliomas, aes(x=X1, y=X2, col=Type)) +
#     geom_point()  +
#     geom_text(data = sDF_gliomas,aes(x=X1,y=X2,label=TypeSurvival), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and selected 1p19q inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()
# 
# message('saving tsneplot_gliomas_ID in ', ptsneplot_gliomas_ID)
# png(file = ptsneplot_gliomas_ID, width = 1440, height = 1280) #width = 1846,height = 991
# ggplot(sDF_gliomas, aes(x=X1, y=X2, col=TypeSurvival)) +
#     geom_point()  +
#     geom_text(data = sDF_gliomas,aes(x=X1,y=X2,label=ID), show.legend = F) +
#     theme_bw() + 
#     theme(legend.position = "right") +
#     guides(col=guide_legend(ncol=1)) +
#     ggtitle(paste0("t-Distributed Stochastic Neighbor Embedding(tSNE) based on the principal component analysis(pca) of all probes of the 1p/19q cohort and selected 1p/19q inhouse samples"),
#             subtitle= paste0('The number of principal components calculated for the pca is ', kpca,'. The tSNE is exact with alpha=0 and the perplexity is set to ', perp,'.')
#     )
# dev.off()





# # # # Multidimensional scaling plot
# 
# 
# par(mfrow= c(1,2))
# m = mdsPlot(betas_gliomas, sampNames = sDF_gliomas[colnames(RGset_gliomas),'TypeSurvival'], numPositions=32000)
# l = mdsPlot(betas_gliomas, sampNames = sDF_gliomas[colnames(RGset_gliomas),'ID'], numPositions=32000)
# 
# 
# par(mfrow= c(1,2))
# m = mdsPlot(RGset_cohort, sampNames = DFclinical_cohort[colnames(RGset_cohort),'TypeSurvival'], numPositions=32000)
# l = mdsPlot(RGset_cohort, sampNames = DFclinical_cohort[colnames(RGset_cohort),'ID'], numPositions=32000)

sink()
