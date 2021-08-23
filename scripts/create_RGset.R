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
# sentrix <- './../../LGG_Methylation/data/utrecht_methy/205061430013_R06C01_Grn.idat'
# sentrix <- list.files('./../../LGG_Methylation/data/allsamples/','*_Grn.idat')

# output file paths
pRGset <- snakemake@output[["RGset"]]
# pRGset <- './../results/RGset/205061430013_R06C01_RGset.RDS'

# parameters
dir_rgset <- snakemake@params[['dir_rgset']]
# dir_rgset <- './../../LGG_Methylation/data/utrecht_methy/'
if(!dir.exists(dir_rgset))(dir.create(dir_rgset))

dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- './../../LGG_Methylation/mnp_training/'

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# read Sentrix files into Red-Green channel
sentrix <- gsub("_Grn.idat","", sentrix) 
RGset <- read.metharray(sentrix, force=TRUE, verbose=TRUE)

message('saving RGset in ', pRGset)
saveRDS(RGset, file = pRGset)

sink()
