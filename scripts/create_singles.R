###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
} else{
    library(minfi)
}


##################
# input file paths
sentrix <- snakemake@input[["sentrix"]]

pRGset <- snakemake@input[["RGset"]]
pMset <- snakemake@input[["Mset"]]
pMset_mnp_filtered <- snakemake@input[["Mset_mnp_filtered"]]
pMset_raw <- snakemake@input[["Mset_raw"]]
pMset_noob <- snakemake@input[["Mset_noob"]]
# sentrix <- './../../LGG_Methylation/data/utrecht_methy/205061430013_R06C01_Grn.idat'

# output file paths
psingleRGset <- snakemake@output[["singleRGset"]]
psingleMset <- snakemake@output[["singleMset"]]
psingleMset_mnp_filtered <- snakemake@output[["singleMset_mnp_filtered"]]
psingleMset_raw <- snakemake@output[["singleMset_raw"]]
psingleMset_noob <- snakemake@output[["singleMset_noob"]]
# pRGset <- './../results/RGset/205061430013_R06C01_RGset.RDS'
# pRGset_mnp_filtered <- './../results/RGset/205061430013_R06C01_RGset_mnp_filtered.RDS'

# parameters
dir_idats <- snakemake@params[['dir_idats']]
dir_rgset <- snakemake@params[['dir_rgset']]
dir_rgset_singles <- snakemake@params[['dir_rgset_singles']]
dir_mset <- snakemake@params[['dir_mset']]
dir_mset_singles <- snakemake@params[['dir_mset_singles']]
# dir_rgset <- './../../LGG_Methylation/data/utrecht_methy/'
if(!dir.exists(dir_rgset))(dir.create(dir_rgset))
if(!dir.exists(dir_rgset_singles))(dir.create(dir_rgset_singles))
if(!dir.exists(dir_mset))(dir.create(dir_mset))
if(!dir.exists(dir_mset_singles))(dir.create(dir_mset_singles))


log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# read input
RGset <- readRDS(pRGset)
Mset <- readRDS(pMset)
Mset_mnp_filtered <- readRDS(pMset_mnp_filtered)
Mset_raw <- readRDS(pMset_raw)
Mset_noob <- readRDS(pMset_noob)

# select sentrix
sentrix <- gsub(dir_idats,"", sentrix)
sentrix <- gsub("_Grn.idat","", sentrix) 

# select sentrix from full RGset
singleRGset <- RGset[,sentrix]
message('saving singleRGset in ', psingleRGset)
saveRDS(singleRGset, file = psingleRGset)

# select sentrix from full Mset
singleMset <- Mset[,sentrix]
message('saving singleMset in ', psingleMset)
saveRDS(singleMset, file = psingleMset)

# select sentrix from full Mset_mnp_filtered
singleMset_mnp_filtered <- Mset_mnp_filtered[,sentrix]
message('saving singleMset_mnp_filtered in ', psingleMset_mnp_filtered)
saveRDS(singleMset_mnp_filtered, file = psingleMset_mnp_filtered)

# select sentrix from full Mset_raw
singleMset_raw <- Mset_raw[,sentrix]
message('saving singleMset_raw in ', psingleMset_raw)
saveRDS(singleMset_raw, file = psingleMset_raw)

# select sentrix from full Mset_noob
singleMset_noob <- Mset_noob[,sentrix]
message('saving singleMset_noob in ', psingleMset_noob)
saveRDS(singleMset_noob, file = psingleMset_noob)

sink()
sink()
