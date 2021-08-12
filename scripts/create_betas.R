###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
    suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
    suppressMessages(library(RSpectra))
} else{
    library(minfi)
    library(IlluminaHumanMethylationEPICmanifest)
    library(RSpectra)
}

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

##################
# input file paths
pMset <- snakemake@input[["Mset"]]
pMset_mnp_filtered <- snakemake@input[["Mset_mnp_filtered"]]
# pMset <- './../results/Mset/205061430013_R06C01_Mset.RDS'
# pMset_mnp_filtered <- './../results/Mset/205061430013_R06C01_Mset_mnp_filtered.RDS'

# output file paths
pbetas <- snakemake@output[["betas"]]
pbetas_mnp_filtered <- snakemake@output[["betas_mnp_filtered"]]
# pbetas_mnp_filtered_sd <- snakemake@output[["betas_mnp_filtered_sd"]]
# pbetas <- './../results/betas/205061430013_R06C01_betas.RDS'
# pbetas_mnp_filtered <- './../results/betas/205061430013_R06C01_betas_mnp_filtered.RDS'
# pbetas_mnp_filtered_sd <- './../results/betas/205061430131_R01C01_betas_mnp_filtered_sd.RDS'

# parameters
# sd_filter <- as.integer(snakemake@params[['sd_filter']])
dir_betas <- snakemake@params[["dir_betas"]]
# sd_filter <- 32000

if(!dir.exists(dir_betas))(dir.create(dir_betas))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# read input
Mset <- readRDS(pMset)
Mset_mnp_filtered <- readRDS(pMset_mnp_filtered)

# calculate betas
betas <- getBeta(Mset, offset = 100)
betas <- as.data.frame(t(betas))

message('saving betas in ', pbetas)
saveRDS(betas, file = pbetas)

betas_mnp_filtered <- getBeta(Mset_mnp_filtered, offset = 100)
betas_mnp_filtered <- as.data.frame(t(betas_mnp_filtered))

message('saving betas_mnp_filtered in ', pbetas_mnp_filtered)
saveRDS(betas_mnp_filtered, file = pbetas_mnp_filtered)

# # order betas on most variance
# message('filter betas on the sd_filter methylation sites with the largest spread')
# betas_mnp_filtered_sd <- betas_mnp_filtered[,order(-apply(betas_mnp_filtered,2,sd))[1:sd_filter]]
# 
# message('saving betas_mnp_filtered_sd in ', pbetas_mnp_filtered_sd)
# saveRDS(betas_mnp_filtered_sd, file = pbetas_mnp_filtered_sd)
# 

sink()
sink()
