
###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    # suppressMessages(library(minfi))
    suppressMessages(library(openxlsx))
    suppressMessages(library(dplyr))
    suppressMessages(library(gplots))
    suppressMessages(library(Rtsne))
    suppressMessages(library(ggplot2))
    suppressMessages(library(purrr))
    suppressMessages(library(RSpectra))
    suppressMessages(library(randomForest))
    suppressMessages(library(pROC))
    suppressMessages(library(RColorBrewer))
} else{
    # library(minfi)
    library(openxlsx)
    library(dplyr)
    library(gplots)
    library(Rtsne)
    library(ggplot2)
    library(purrr)
    library(RSpectra)
    library(randomForest)
    library(pROC)
    library(RColorBrewer)
}

##################
# # # input file paths
pbetas <- snakemake@input[["betas_cohort"]]
pbetas_gliomas <- snakemake@input[["betas_gliomas"]]
pbetas_inhouse <- snakemake@input[["betas_inhouse"]]
pDFclinical_cohort <-        snakemake@input[["DFclinical_cohort"]]
pDFclinical_gliomas <-      snakemake@input[["DFclinical_gliomas"]]
pDFclinical_inhouse <-      snakemake@input[["DFclinical_inhouse"]]
# pbetas <- './../results/betas/betas.RDS'
# # # output file paths
# # betas
pBetas_gCIMP <- snakemake@output[["betas_gCIMP"]]
pBetas_gCIMP_short <-       snakemake@output[["betas_gCIMP_short"]]
pBetas_gCIMP_long <-        snakemake@output[["betas_gCIMP_long"]]
# pBetas_gCIMP_thresholded <- snakemake@output[["betas_gCIMP_thresholded"]]
# pBetas_gCIMP_short_thresholded <- snakemake@output[["betas_gCIMP_short_thresholded"]]
# pBetas_gCIMP_long_thresholded <- snakemake@output[["betas_gCIMP_long_thresholded"]]
pBetas_glass_treatment_related_620probes <-         snakemake@output[["betas_glass_treatment_related_620probes"]]
pBetas_glass_treatment_related_620probes_short <-   snakemake@output[["betas_glass_treatment_related_620probes_short"]]
pBetas_glass_treatment_related_620probes_long <-    snakemake@output[["betas_glass_treatment_related_620probes_long"]]
pBetas_glass_hypermodulator_342probes <-            snakemake@output[["betas_glass_hypermodulator_342probes"]]
pBetas_glass_hypermodulator_342probes_short <-      snakemake@output[["betas_glass_hypermodulator_342probes_short"]]
pBetas_glass_hypermodulator_342probes_long <-       snakemake@output[["betas_glass_hypermodulator_342probes_long"]]
#     pBetas_gCIMP <-                     './../results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS'
#     pBetas_gCIMP_short <-               './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_short.RDS'
#     pBetas_gCIMP_long <-                './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_long.RDS'
# pBetas_gCIMP_thresholded <-         './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_thresholded.RDS'
#     pBetas_gCIMP_short_thresholded <-   './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_short_thresholded.RDS'
#     pBetas_gCIMP_long_thresholded <-    './../results/analysis_probes/gCIMP_7probes/betas_gCIMP_long_thresholded.RDS'

# # heatmap
pheatmap_gCIMP_cohort <-                snakemake@output[['heatmap_gCIMP_cohort']]
pheatmap_gCIMP_cohort_short <-          snakemake@output[['heatmap_gCIMP_cohort_short']]
pheatmap_gCIMP_cohort_long <-           snakemake@output[['heatmap_gCIMP_cohort_long']]
# pheatmap_gCIMP_cohort_thresholded <-    snakemake@output[['heatmap_gCIMP_cohort_thresholded']]
# pheatmap_gCIMP_cohort_thresholded_short <-  snakemake@output[['heatmap_gCIMP_cohort_thresholded_short.png']]
# pheatmap_gCIMP_cohort_thresholded_long <-   snakemake@output[['heatmap_gCIMP_cohort_thresholded_long.png']]
#     pheatmap_gCIMP_cohort <-        paste0(dir_gCIMP,'heatmap_gCIMP_cohort.png')
#     pheatmap_gCIMP_cohort_short <-  paste0(dir_gCIMP,'heatmap_gCIMP_cohort_short.png')
#     pheatmap_gCIMP_cohort_long <-   paste0(dir_gCIMP,'heatmap_gCIMP_cohort_long.png')
#     pheatmap_gCIMP_cohort_thresholded <-        paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded.png')
#     pheatmap_gCIMP_cohort_thresholded_short <-  paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded_short.png')
#     pheatmap_gCIMP_cohort_thresholded_long <-   paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded_long.png')
pheatmap_treatment_related_cohort <-        snakemake@output[['heatmap_treatment_related_cohort']]
pheatmap_treatment_related_cohort_short <-  snakemake@output[['heatmap_treatment_related_cohort_short']]
pheatmap_treatment_related_cohort_long <-   snakemake@output[['heatmap_treatment_related_cohort_long']]

pheatmap_hypermodulator_cohort <-        snakemake@output[['heatmap_hypermodulator_cohort']]
pheatmap_hypermodulator_cohort_short <-  snakemake@output[['heatmap_hypermodulator_cohort_short']]
pheatmap_hypermodulator_cohort_long <-   snakemake@output[['heatmap_hypermodulator_cohort_long']]
# # heatmaps RDS
pheatmapRDS_gCIMP_cohort <-                snakemake@output[['heatmapRDS_gCIMP_cohort']]
pheatmapRDS_gCIMP_cohort_short <-          snakemake@output[['heatmapRDS_gCIMP_cohort_short']]
pheatmapRDS_gCIMP_cohort_long <-           snakemake@output[['heatmapRDS_gCIMP_cohort_long']]
# pheatmapRDS_gCIMP_cohort_thresholded <-    snakemake@output[['heatmapRDS_gCIMP_cohort_thresholded']]
# pheatmap_gCIMP_cohort_thresholded_short <-  snakemake@output[['heatmap_gCIMP_cohort_thresholded_short.png']]
# pheatmap_gCIMP_cohort_thresholded_long <-   snakemake@output[['heatmap_gCIMP_cohort_thresholded_long.png']]
#     pheatmap_gCIMP_cohort <-        paste0(dir_gCIMP,'heatmap_gCIMP_cohort.png')
#     pheatmap_gCIMP_cohort_short <-  paste0(dir_gCIMP,'heatmap_gCIMP_cohort_short.png')
#     pheatmap_gCIMP_cohort_long <-   paste0(dir_gCIMP,'heatmap_gCIMP_cohort_long.png')
#     pheatmap_gCIMP_cohort_thresholded <-        paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded.png')
#     pheatmap_gCIMP_cohort_thresholded_short <-  paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded_short.png')
#     pheatmap_gCIMP_cohort_thresholded_long <-   paste0(dir_gCIMP,'heatmap_gCIMP_cohort_thresholded_long.png')
pheatmapRDS_treatment_related_cohort <-        snakemake@output[['heatmapRDS_treatment_related_cohort']]
pheatmapRDS_treatment_related_cohort_short <-  snakemake@output[['heatmapRDS_treatment_related_cohort_short']]
pheatmapRDS_treatment_related_cohort_long <-   snakemake@output[['heatmapRDS_treatment_related_cohort_long']]

pheatmapRDS_hypermodulator_cohort <-        snakemake@output[['heatmapRDS_hypermodulator_cohort']]
pheatmapRDS_hypermodulator_cohort_short <-  snakemake@output[['heatmapRDS_hypermodulator_cohort_short']]
pheatmapRDS_hypermodulator_cohort_long <-   snakemake@output[['heatmapRDS_hypermodulator_cohort_long']]

# # tSNE
# ptsne_cohort_gliomas_RDS <- snakemake@output[['tsne_cohort_gCIMP_gliomas_RDS']]
# ptsne_cohort_gliomas_Type <- snakemake@output[['tsne_cohort_gCIMP_gliomas_Type']]
# ptsne_cohort_gliomas_TypeSurvival <- snakemake@output[['tsne_cohort_gCIMP_gliomas_TypeSurvival']]
#     ptsne_gCIMP_gliomas_RDS <-              paste0(dir_gCIMP, 'tsne_gCIMP_gliomas.RDS')
#     ptsne_gCIMP_gliomas_Type <-             paste0(dir_gCIMP, 'tsne_gCIMP_gliomas_Type.png')
#     ptsne_gCIMP_gliomas_TypeSurvival <-     paste0(dir_gCIMP, 'tsne_gCIMP_gliomas_TypeSurvival.png')

# # # parameters
dir_analysis_probes <- snakemake@params[['dir_analysis_probes']]
dir_gCIMP <- snakemake@params[['dir_gCIMP']]
dir_glass <- snakemake@params[['dir_glass']]
dir_glass_treatment_related_probes <-   snakemake@params[['dir_glass_treatment_related_probes']]
dir_glass_hypermodulator_probes <-      snakemake@params[['dir_glass_hypermodulator_probes']]
pprobeset_gCIMP <- snakemake@params[['probeset_gCIMP']]
pprobeset_GLASS <- snakemake@params[['probeset_GLASS']]
probeset_GLASS_treatment_related_sheetname <-   snakemake@params[['probeset_GLASS_treatment_related_sheetname']]
probeset_GLASS_hypermodulator_sheetname <-      snakemake@params[['probeset_GLASS_hypermodulator_sheetname']]
#     dir_analysis_probes <- './../results/analysis_probes/' 
#     dir_gCIMP <-                            paste0(dir_analysis_probes,'gCIMP_7probes/')
#     dir_glass <-                            paste0(dir_analysis_probes,'GLASS/')
#     dir_glass_treatment_related_probes <-   paste0(dir_glass,'treatment_related_620probes/')
#     dir_glass_hypermodulator_probes <-      paste0(dir_glass,'hypermodulator_342probes/')
#     pprobeset_gCIMP <- './../../LGG_Methylation/ProbeSet-gCIMP.xlsx'
#     pprobeset_GLASS <- './../../LGG_Methylation/ProbeSet-GLASS.xlsx'
#     probeset_GLASS_treatment_related_sheetname <-   'Treatment-related (N=620)'
#     probeset_GLASS_hypermodulator_sheetname <-      'IDHmut Hypermutator (N=342)'

if(!dir.exists(dir_analysis_probes))(dir.create(dir_analysis_probes))
if(!dir.exists(dir_gCIMP))(dir.create(dir_gCIMP))
if(!dir.exists(dir_glass))(dir.create(dir_glass))
if(!dir.exists(dir_glass_treatment_related_probes))(dir.create(dir_glass_treatment_related_probes))
if(!dir.exists(dir_glass_hypermodulator_probes))(dir.create(dir_glass_hypermodulator_probes))

pinhouse <- snakemake@params[['inhouse']]
pcohort <- snakemake@params[['overview']]
# pinhouse <- './../../LGG_Methylation/Yongsoo VUMC methylation arrat samples sep 2020.xlsx'
# pcohort <- './../../LGG_Methylation/Overview_oligodendroglioma_samples_long_vs_short_PFS_2021-06-10.xlsx'


dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# # # Read input
betas <- readRDS(pbetas)
betas_gliomas <- readRDS(pbetas_gliomas)
betas_inhouse <- readRDS(pbetas_inhouse)
DFclinical_cohort <- readRDS(pDFclinical_cohort)
DFclinical_gliomas <- readRDS(pDFclinical_gliomas)
DFclinical_inhouse <- readRDS(pDFclinical_inhouse)


# # # select probes gCIMP
# cgs1 <- ('cg09732711', 'cg09326832', 'cg24665265', 'cg06220958', 'cg10245915', 'cg11689625', 'cg11799650')
cgs1_ordered <- c('cg09732711', 'cg24665265', 'cg10245915', 'cg11799650', 'cg09326832', 'cg06220958', 'cg11689625')
cgs1 <- read.xlsx(pprobeset_gCIMP)
cgs1 <- cgs1[,1]

# # # define betas for full cohort and gCIMP probes
betas_probes7_cohort <- betas[,cgs1]
message('saving betas for full cohort and gCIMP probes in ', pBetas_gCIMP)
saveRDS(betas_probes7_cohort, file = pBetas_gCIMP)

# # # define betas for short cohort and gCIMP probes
betas_probes7_cohort_short <- betas_probes7_cohort[DFclinical_cohort$TypeSurvival=='short',]
message('saving betas for short cohort and gCIMP probes in ', pBetas_gCIMP_short)
saveRDS(betas_probes7_cohort_short, file = pBetas_gCIMP_short)

# # # define betas for long cohort and gCIMP probes
betas_probes7_cohort_long <- betas_probes7_cohort[DFclinical_cohort$TypeSurvival=='long',]
message('saving betas for long cohort and gCIMP probes in ', pBetas_gCIMP_long)
saveRDS(betas_probes7_cohort_long, file = pBetas_gCIMP_long)

# # # cut-off values for betas for gCIMP probes
# cg09732711 (0.7), cg09326832 (0.28), cg24665265 (0.67),cg06220958 (0.17), cg10245915 (0.12), cg11689625 (0.31), and cg11799650 (0.49)

# cgs1 <- c('cg09732711', 'cg09326832', 'cg24665265', 'cg06220958', 'cg10245915', 'cg11689625', 'cg11799650')
cutoffs <- c(0.7,0.28,0.67,0.17,0.12,0.31,0.49)
betas_probes7_cohort_thresholded <- betas_probes7_cohort
for (i in seq_along(cutoffs)){
    thres <- cutoffs[i]
    probe <- cgs1[i]
    values <- betas_probes7_cohort[,probe] >= thres
    betas_probes7_cohort_thresholded[,probe] <- as.numeric(values)
}
# message('saving thresholded betas for full cohort of gCIMP probes in ', pBetas_gCIMP_thresholded)
# saveRDS(betas_probes7_cohort_thresholded, file = pBetas_gCIMP_thresholded)


# # # functions for heatmaps
# betas_probes7_cohort = readRDS('./../results/analysis_probes/gCIMP_7probes/betas_gCIMP.RDS')
# DFclinical_cohort = readRDS('./../../output-lgg/results/DFclinical_cohort.RDS')
# betas_probes7_cohort = betas_probes7_cohort[setdiff(rownames(betas_probes7_cohort),'205061430022_R06C01'),]
fdist <- function(x) as.dist(1-cor(t(x)))
fclus <- function(x) hclust(x,method= "ward.D2")

# d <- fdist(betas_probes7_cohort)
# fit <- fclus(d)

colAnn <- colorpanel(2,low="red",high="green")
# DFclinical_cohort = DFclinical_cohort[setdiff(rownames(DFclinical_cohort),'205061430022_R06C01'),]
DFclinical_cohort$colAnn <- ifelse(DFclinical_cohort$TypeSurvival == 'short', colAnn[1], colAnn[2])

# # make heatmap for betas
# order heatmap as in gCIMP figure 5A
message('saving heatmap for full cohort and gCIMP probes in ', pheatmap_gCIMP_cohort)
png(file = pheatmap_gCIMP_cohort, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
# png(width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p1 = heatmap.2(as.matrix(betas_probes7_cohort[,cgs1_ordered]),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort$colAnn,
          labRow = DFclinical_cohort[rownames(betas_probes7_cohort[,cgs1_ordered]),'ID'],
          trace = 'none',
          Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes7_cohort[,cgs1_ordered]),'ID'],DFclinical_cohort[rownames(betas_probes7_cohort[,cgs1_ordered]),'TypeSurvival']),
          main=paste('heatmap of 7 beta values of methylation probes\n for full cohort of',length(rownames(betas_probes7_cohort[,cgs1_ordered])),'samples of 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
legend("topright",      
    legend = unique(DFclinical_cohort$TypeSurvival),
    title = 'TypeSurvival',
    col = unique(DFclinical_cohort$colAnn), 
    lty= 1,             
    lwd = 5,           
    cex=1
    )
print(p1)
dev.off()
message('saving heatmap RDS for full cohort and gCIMP probes in ', pheatmapRDS_gCIMP_cohort)
saveRDS(p1, file = pheatmapRDS_gCIMP_cohort)


message('saving heatmap for short survivors of cohort and gCIMP probes in ', pheatmap_gCIMP_cohort_short)
png(file = pheatmap_gCIMP_cohort_short, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p2 = heatmap.2(as.matrix(betas_probes7_cohort_short[,cgs1_ordered]),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          #           RowSideColors = DFclinical_cohort[rownames(betas_probes7_cohort_short),'colAnn'],
          trace = 'none',
          Colv = FALSE,
          srtCol = 45,
          labRow = DFclinical_cohort[rownames(betas_probes7_cohort_short[,cgs1_ordered]),'ID'],
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes7_cohort_short[,cgs1_ordered]),'ID'],DFclinical_cohort[rownames(betas_probes7_cohort_short[,cgs1_ordered]),'TypeSurvival']),
          main=paste('heatmap of beta values of 7 methylation probes\n for full cohort of',length(rownames(betas_probes7_cohort_short)),'samples of short survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p2)
dev.off()
message('saving heatmap RDS for short survivors of cohort and gCIMP probes in ', pheatmapRDS_gCIMP_cohort_short)
saveRDS(p2, file = pheatmapRDS_gCIMP_cohort_short)

message('saving heatmap for long survivors of cohort and gCIMP probes in ', pheatmap_gCIMP_cohort_long)
png(file = pheatmap_gCIMP_cohort_long, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p3 = heatmap.2(as.matrix(betas_probes7_cohort_long[,cgs1_ordered]),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          trace = 'none',
          Colv = FALSE,
          srtCol = 45,
          labRow = DFclinical_cohort[rownames(betas_probes7_cohort_long[,cgs1_ordered]),'ID'],
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes7_cohort_long[,cgs1_ordered]),'ID'],DFclinical_cohort[rownames(betas_probes7_cohort_long[,cgs1_ordered]),'TypeSurvival']),
          main=paste('heatmap of beta values of 7 methylation probes\n for full cohort of',length(rownames(betas_probes7_cohort_long)),'samples of long survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p3)
dev.off()
message('saving heatmap RDS for long survivors of cohort and gCIMP probes in ', pheatmapRDS_gCIMP_cohort_long)
saveRDS(p3, file = pheatmapRDS_gCIMP_cohort_long)

# # # # heatmaps for thresholded gCIMP betas
# message('saving heatmap for full cohort and thresholded gCIMP probes in ', pheatmap_gCIMP_cohort_thresholded)
# png(file = pheatmap_gCIMP_cohort_thresholded, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
# p4 = heatmap.2(as.matrix(betas_probes7_cohort_thresholded[,cgs1_ordered]),
#           col = rev(brewer.pal(n=11,'RdYlBu')),
#           hclustfun = fclus, distfun = fdist,
#           RowSideColors = DFclinical_cohort$colAnn,
#           labRow = DFclinical_cohort[rownames(betas_probes7_cohort_thresholded[,cgs1_ordered]),'ID'],
#           trace = 'none',
#           Colv = FALSE,
#           srtCol = 45,
#           #           labRow = paste(DFclinical_cohort[rownames(betas_probes7_cohort_thresholded[,cgs1_ordered]),'ID'],DFclinical_cohort[rownames(betas_probes7_cohort_thresholded[,cgs1_ordered]),'TypeSurvival']),
#           main=paste('heatmap of 7 thresholded beta values of methylation probes\n for full cohort of',length(rownames(betas_probes7_cohort_thresholded[,cgs1_ordered])),'samples of 1p/19q oligodendroglioma'),
#           margin=c(10,20)
# )
# legend("topright",      
#     legend = unique(DFclinical_cohort$TypeSurvival),
#     title = 'TypeSurvival',
#     col = unique(DFclinical_cohort$colAnn), 
#     lty= 1,             
#     lwd = 5,           
#     cex=.7
#     )
# print(p4)
# dev.off()
# message('saving heatmap RDS for full cohort and thresholded gCIMP probes in ', pheatmapRDS_gCIMP_cohort_thresholded)
# saveRDS(p4, file = pheatmapRDS_gCIMP_cohort_thresholded)



# # # select probes glass treatment_related_620probes
cgs1 <- read.xlsx(pprobeset_GLASS, sheet=probeset_GLASS_treatment_related_sheetname)
cgs1 <- cgs1[,1]
betas_probes620 <- betas[,cgs1]
DFclinical_cohort <- DFclinical_cohort
betas_probes620_cohort <- betas_probes620[rownames(DFclinical_cohort),]

# # # define betas for full cohort and treatment related probes
betas_probes620_cohort <- betas_probes620_cohort
message('saving betas for full cohort and treatment related probes in ', pBetas_glass_treatment_related_620probes)
saveRDS(betas_probes620_cohort, file = pBetas_glass_treatment_related_620probes)

# # # define betas for short cohort and treatment related probes
betas_probes620_cohort_short <- betas_probes620_cohort[DFclinical_cohort$TypeSurvival=='short',]
message('saving betas for short cohort and treatment related probes in ', pBetas_glass_treatment_related_620probes_short)
saveRDS(betas_probes620_cohort_short, file = pBetas_glass_treatment_related_620probes_short)

# # # define betas for long cohort and treatment related probes
betas_probes620_cohort_long <- betas_probes620_cohort[DFclinical_cohort$TypeSurvival=='long',]
message('saving betas for long cohort and treatment related probes in ', pBetas_glass_treatment_related_620probes_long)
saveRDS(betas_probes620_cohort_long, file = pBetas_glass_treatment_related_620probes_long)

# # make heatmap for betas
# full
# rownames(betas_probes620_cohort) <- paste(DFclinical_cohort[rownames(betas_probes620_cohort),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort),'TypeSurvival'], 'survivor')
message('saving heatmap for full cohort and treatment related probes in ', pheatmap_treatment_related_cohort)
png(file = pheatmap_treatment_related_cohort, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p1 = heatmap.2(as.matrix(betas_probes620_cohort),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort$colAnn,
          labRow = DFclinical_cohort[rownames(betas_probes620_cohort),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes620_cohort),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort),'TypeSurvival']),
          main=paste('heatmap of 620 beta values of methylation probes\n for full cohort of',length(rownames(betas_probes620_cohort)),'samples of 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
legend("topright",      
    legend = unique(DFclinical_cohort$TypeSurvival),
    title = 'TypeSurvival',
    col = unique(DFclinical_cohort$colAnn), 
    lty= 1,             
    lwd = 5,           
    cex=1
    )
print(p1)
dev.off()
message('saving heatmap RDS for full cohort and treatment related probes in ', pheatmapRDS_treatment_related_cohort)
saveRDS(p1, file = pheatmapRDS_treatment_related_cohort)
# short
rownames(betas_probes620_cohort_short) <- paste(DFclinical_cohort[rownames(betas_probes620_cohort_short),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort_short),'TypeSurvival'], 'survivor')
message('saving heatmap for short survivors of cohort and treatment related probes in ', pheatmap_treatment_related_cohort_short)
png(file = pheatmap_treatment_related_cohort_short, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p2 = heatmap.2(as.matrix(betas_probes620_cohort_short),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort[rownames(betas_probes620_cohort_short),'colAnn'],
          labRow = DFclinical_cohort[rownames(betas_probes620_cohort_short),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes620_cohort_short),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort_short),'TypeSurvival']),
          main=paste('heatmap of beta values of 620 methylation probes\n for cohort of',length(rownames(betas_probes620_cohort_short)),'samples of short survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p2)
dev.off()
message('saving heatmap RDS for short survivors of cohort and treatment related probes in ', pheatmapRDS_treatment_related_cohort_short)
saveRDS(p2, file = pheatmapRDS_treatment_related_cohort_short)
# long
rownames(betas_probes620_cohort_long) <- paste(DFclinical_cohort[rownames(betas_probes620_cohort_long),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort_long),'TypeSurvival'], 'survivor')
message('saving heatmap for long survivors of cohort and treatment related probes in ', pheatmap_treatment_related_cohort_long)
png(file = pheatmap_treatment_related_cohort_long, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p3 = heatmap.2(as.matrix(betas_probes620_cohort_long),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort[rownames(betas_probes620_cohort_long),'colAnn'],
          labRow = DFclinical_cohort[rownames(betas_probes620_cohort_long),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes620_cohort_long),'ID'],DFclinical_cohort[rownames(betas_probes620_cohort_long),'TypeSurvival']),
          main=paste('heatmap of beta values of 620 methylation probes\n for cohort of',length(rownames(betas_probes620_cohort_long)),'samples of long survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p3)
dev.off()
message('saving heatmap RDS for long survivors of cohort and treatment related probes in ', pheatmapRDS_treatment_related_cohort_long)
saveRDS(p3, file = pheatmapRDS_treatment_related_cohort_long)


# # # select probes glass hypermodulator_342probes
cgs1 <- read.xlsx(pprobeset_GLASS, sheet=probeset_GLASS_hypermodulator_sheetname)
cgs1 <- cgs1[,1]
betas_probes342 <- betas[,cgs1]
DFclinical_cohort <- DFclinical_cohort
betas_probes342_cohort <- betas_probes342[rownames(DFclinical_cohort),]

# # # define betas for full cohort and hypermodulator probes
betas_probes342_cohort <- betas_probes342_cohort
message('saving betas for full cohort and hypermodulator probes in ', pBetas_glass_hypermodulator_342probes)
saveRDS(betas_probes342_cohort, file = pBetas_glass_hypermodulator_342probes)

# # # define betas for short cohort and hypermodulator probes
betas_probes342_cohort_short <- betas_probes342_cohort[DFclinical_cohort$TypeSurvival=='short',]
message('saving betas for short cohort and hypermodulator probes in ', pBetas_glass_hypermodulator_342probes_short)
saveRDS(betas_probes342_cohort_short, file = pBetas_glass_hypermodulator_342probes_short)

# # # define betas for long cohort and hypermodulator probes
betas_probes342_cohort_long <- betas_probes342_cohort[DFclinical_cohort$TypeSurvival=='long',]
message('saving betas for long cohort and hypermodulator probes in ', pBetas_glass_hypermodulator_342probes_long)
saveRDS(betas_probes342_cohort_long, file = pBetas_glass_hypermodulator_342probes_long)

# # make heatmap for betas
# full
# rownames(betas_probes342_cohort) <- paste(DFclinical_cohort[rownames(betas_probes342_cohort),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort),'TypeSurvival'], 'survivor')
message('saving heatmap for full cohort and hypermodulator probes in ', pheatmap_hypermodulator_cohort)
png(file = pheatmap_hypermodulator_cohort, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p1 = heatmap.2(as.matrix(betas_probes342_cohort),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort$colAnn,
          labRow = DFclinical_cohort[rownames(betas_probes342_cohort),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes342_cohort),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort),'TypeSurvival']),
          main=paste('heatmap of 342 beta values of methylation probes\n for full cohort of',length(rownames(betas_probes342_cohort)),'samples of 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
legend("topright",      
    legend = unique(DFclinical_cohort$TypeSurvival),
    col = unique(DFclinical_cohort$colAnn), 
    lty= 1,             
    lwd = 5,           
    cex=.7
    )
print(p1)
dev.off()
message('saving heatmap RDS for full cohort and hypermodulator probes in ', pheatmapRDS_hypermodulator_cohort)
saveRDS(p1, file = pheatmapRDS_hypermodulator_cohort)
# short
rownames(betas_probes342_cohort_short) <- paste(DFclinical_cohort[rownames(betas_probes342_cohort_short),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort_short),'TypeSurvival'], 'survivor')
message('saving heatmap for short survivors of cohort and hypermodulator probes in ', pheatmap_hypermodulator_cohort_short)
png(file = pheatmap_hypermodulator_cohort_short, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p2 = heatmap.2(as.matrix(betas_probes342_cohort_short),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort[rownames(betas_probes620_cohort_short),'colAnn'],
          labRow = DFclinical_cohort[rownames(betas_probes620_cohort_short),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes342_cohort_short),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort_short),'TypeSurvival']),
          main=paste('heatmap of beta values of 342 methylation probes\n for full cohort of',length(rownames(betas_probes342_cohort_short)),'samples of short survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p2)
dev.off()
message('saving heatmap RDS for short survivors of cohort and hypermodulator probes in ', pheatmapRDS_hypermodulator_cohort_short)
saveRDS(p2, file = pheatmapRDS_hypermodulator_cohort_short)
# long
rownames(betas_probes342_cohort_long) <- paste(DFclinical_cohort[rownames(betas_probes342_cohort_long),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort_long),'TypeSurvival'], 'survivor')
message('saving heatmap for long survivors of cohort and hypermodulator probes in ', pheatmap_hypermodulator_cohort_long)
png(file = pheatmap_hypermodulator_cohort_long, width = 1000, height = 1000, pointsize = 14) #width = 1846,height = 991
p3 = heatmap.2(as.matrix(betas_probes342_cohort_long),
          col = rev(brewer.pal(n=11,'RdYlBu')),
          hclustfun = fclus, distfun = fdist,
          RowSideColors = DFclinical_cohort[rownames(betas_probes342_cohort_long),'colAnn'],
          labRow = DFclinical_cohort[rownames(betas_probes342_cohort_long),'ID'],
          trace = 'none',
          #           Colv = FALSE,
          srtCol = 45,
          #           labRow = paste(DFclinical_cohort[rownames(betas_probes342_cohort_long),'ID'],DFclinical_cohort[rownames(betas_probes342_cohort_long),'TypeSurvival']),
          main=paste('heatmap of beta values of 342 methylation probes\n for full cohort of',length(rownames(betas_probes342_cohort_long)),'samples of long survivors with 1p/19q oligodendroglioma'),
          margin=c(10,20)
)
print(p3)
dev.off()
message('saving heatmap RDS for long survivors of cohort and hypermodulator probes in ', pheatmapRDS_hypermodulator_cohort_long)
saveRDS(p3, file = pheatmapRDS_hypermodulator_cohort_long)

sink()
