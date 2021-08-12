###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
    suppressMessages(library(RSpectra))
} else{
    library(minfi)
    library(RSpectra)
}

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- "./../../LGG_Methylation/mnp_training/"

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

##################
# input file paths
pRGset <- snakemake@input[["RGset"]]
# pRGset <- './../results/RGset/205061430013_R06C01_RGset.RDS'

# output file paths
pMset <- snakemake@output[["Mset"]]
pMset_mnp_filtered <- snakemake@output[["Mset_mnp_filtered"]]
pMset_raw <- snakemake@output[["Mset_raw"]]
pMset_noob <- snakemake@output[["Mset_noob"]]
# pMset <- './../results/Mset/205061430013_R06C01_Mset.RDS'
# pMset_mnp_filtered <- './../results/Mset/205061430013_R06C01_Mset_mnp_filtered.RDS'
# pMset_raw <- './../results/Mset/205061430013_R06C01_Mset_raw.RDS'

# parameters
dir_mset <- snakemake@params[["dir_mset"]]
if(!dir.exists(dir_mset))(dir.create(dir_mset))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# read input
RGset <- readRDS(pRGset)

# MNP Illumina normalization
message("running Illumina normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(RGset)

message('saving Mset of inhouse and cohort samples with quality control  in ', pMset)
saveRDS(Mset, file = pMset)

# Raw normalization for CNV
message("running raw normalization ...",Sys.time())
Mset_raw <- preprocessRaw(RGset)

message('saving Mset_raw in ', pMset_raw)
saveRDS(Mset_raw, file = pMset_raw)

# Noob normalization for probeset regression model
message("running noob normalization ...",Sys.time())
Mset_noob <- preprocessNoob(RGset, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference")

message('saving Mset_noob in ', pMset_noob)
saveRDS(Mset_noob, file = pMset_noob)


# probe mnp filtering
message("probe mnp filtering ...",Sys.time())
amb.filter <- read.table(file.path(dir_mnp, "filter","amb_3965probes.vh20151030.txt"),header=F)
epic.filter <- read.table(file.path(dir_mnp, "filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path(dir_mnp, "filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path(dir_mnp, "filter","xy_11551probes.vh20151030.txt"),header=F)
rs.filter <- grep("rs",rownames(Mset))
ch.filter <- grep("ch",rownames(Mset))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(epic.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))
remove <- setdiff(remove, NA)

Mset_mnp_filtered <- Mset[-remove,]

message('saving Mset_mnp_filtered in ', pMset_mnp_filtered)
saveRDS(Mset_mnp_filtered, file = pMset_mnp_filtered)


sink()
sink()
