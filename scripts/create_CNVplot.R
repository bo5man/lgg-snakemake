###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
    suppressMessages(library(conumee))
    suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
    suppressMessages(library(RSpectra))
} else{
    library(minfi)
    library(conumee)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(RSpectra)
}


##################
# input file paths
sentrix <- snakemake@input[["sentrix"]]
pMset_raw <- snakemake@input[["Mset_raw"]]
pDFclinical_cohort <-       snakemake@input[["DFclinical_cohort"]]
pDFclinical_full_cohort <-  snakemake@input[["DFclinical_full_cohort"]]
pDFclinical_gliomas <-      snakemake@input[["DFclinical_gliomas"]]
pDFclinical_inhouse <-      snakemake@input[["DFclinical_inhouse"]]
# sentrix <- '205061430131_R02C01_Grn.idat'
# pMset_raw <- './../results/Mset/Mset_raw.RDS'

# output file paths
pcnv <- snakemake@output[["cnv"]]
pcnvplot <- snakemake@output[["cnvplot"]]
# pcnv <- './../results/cnv/RDS/205061430088_R07C01_cnv.RDS'
# pcnvplot <- './../results/cnv/205061430088_R07C01_cnv.png'

# parmeters
dir_full_cohort <- snakemake@params[['dir_full_cohort']]
dir_mnp <- snakemake@params[['dir_mnp']]
# dir_mnp <- './../../LGG_Methylation/mnp_training/'

dir_cnv <- snakemake@params[["dir_cnv"]]
dir_cnv_rds <- snakemake@params[["dir_cnv_rds"]]
# dir_cnv <- './../results/cnv/'
if(!dir.exists(dir_cnv))(dir.create(dir_cnv))
if(!dir.exists(dir_cnv_rds))(dir.create(dir_cnv_rds))

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

# load conumee annotation object
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_ANNO.vh20150715.RData"))
# load conumee reference male
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_REF-M.vh20150715.RData"))
# load conumee reference female
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_REF-F.vh20150715.RData"))

# read clinical data files
DFclinical_inhouse <- readRDS(pDFclinical_inhouse)
DFclinical_cohort <- readRDS(pDFclinical_cohort)
DFclinical_full_cohort <- readRDS(pDFclinical_full_cohort)
# DFclinical_inhouse = readRDS('./../results/DFclinical_inhouse.RDS')
# DFclinical_cohort = readRDS('./../results/DFclinical_cohort.RDS')

# read raw Mset: no normalization is performed before CNV analysis 
Mset_raw <- readRDS(pMset_raw)

# select current sentrix from Mset
sentrix <- gsub(dir_full_cohort, '', sentrix)
sentrix <- gsub('_Grn.idat', '', sentrix)
sMset <- Mset_raw[,sentrix]

# select probes for EPICarray
# from https://support.bioconductor.org/p/9136138/#9136840
anno@probes <- anno@probes[names(anno@probes) %in% names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]

# get CN data and perform conumee analysis
cndata <- CNV.load(sMset)

# get sex of sample patient
Sex <- DFclinical_full_cohort$SexType[DFclinical_full_cohort$Sentrix_ID == sentrix]
# Sex <- DFclinical_cohort$SexType[DFclinical_cohort$Sentrix_ID == sentrix]

if (Sex == 'F'){
    ref.data = refF.data
}else if (Sex == 'M'){
    ref.data = refM.data
}else{
    stop('stopped in rule CNVplot at sentrix ',sentrix)
}
x <- CNV.fit(cndata, ref.data, anno)
x <- CNV.bin(x)
# x <- CNV.detail(x)
x <- CNV.segment(x)
cnv <- x

message('saving cnv in ', pcnv)
saveRDS(cnv, file = pcnv)

Sentrix_ID <- sentrix
# FullSample <- DFclinical_full_cohort$FullSample[DFclinical_full_cohort$Sentrix_ID == sentrix]
LSS <- DFclinical_full_cohort$LSS[DFclinical_full_cohort$Sentrix_ID == sentrix]
TypeSurvival <- DFclinical_full_cohort$TypeSurvival[DFclinical_full_cohort$Sentrix_ID == sentrix]

# # from QDNAseq
# png(filename = output.file, width = 1846,height = 991)
# plot(data[,i], ylim = c(-1.5,1.5), pointcex = 1.5, pointpch = 16)
# # # plot
  # png(file = pcnvplot, width = 1753, height = 1240) #width = 1846,height = 991
png(file = pcnvplot, width = 1846, height = 991) #width = 1846,height = 991
CNV.genomeplot(cnv, chrY=FALSE, chrX=FALSE, ylim=c(-1.5,1.5),
               #                pointcex=1.5, pointpch=16,
               main = paste('Sentrix:', Sentrix_ID,
                            #                             ', Full Sample Name:',FullSample,
                            ', LSS_ID:', LSS,
                            ', TypeSurvival:', TypeSurvival
               )
)
dev.off()
