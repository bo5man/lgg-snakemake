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
# sentrix <- '205061430055_R04C01_Grn.idat' # E17-022910:LSS22 Female
# pMset_raw <- './../../output-lgg/results/Mset/Mset_raw.RDS'
sentrix <- snakemake@input[["sentrix"]]
pMset_raw <- snakemake@input[["Mset_raw"]]
pDFclinical_full_cohort <-  snakemake@input[["DFclinical_full_cohort"]]

# output file paths
# pcnv <- './../results/cnv/RDS/205061430170_R06C01_cnv.RDS'
# pcnvplot <- './../results/cnv/205061430170_R06C01_cnv.png'
pcnv <- snakemake@output[["cnv"]]
pcnvplot <- snakemake@output[["cnvplot"]]

# parmeters
dir_mnp <- './../../LGG_Methylation/mnp_training/'
# dir_full_cohort <- './../../LGG_Methylation/data/utrecht_methy/'
dir_full_cohort <- snakemake@params[['dir_full_cohort']]
dir_mnp <- snakemake@params[['dir_mnp']]

dir_cnv <- snakemake@params[["dir_cnv"]]
dir_cnv_rds <- snakemake@params[["dir_cnv_rds"]]
# dir_cnv <- './../results/cnv/'
if(!dir.exists(dir_cnv))(dir.create(dir_cnv))
if(!dir.exists(dir_cnv_rds))(dir.create(dir_cnv_rds))

source(file.path(dir_mnp,"R","MNPprocessIDAT_functions.R"))
source(file.path(dir_mnp,"R","RSpectra_pca.R"))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE, type='output')
sink(log, append=T, split=FALSE, type='message')
##################

# load conumee annotation object
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_ANNO.vh20150715.RData"))
# load conumee reference male
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_REF-M.vh20150715.RData"))
# load conumee reference female
load(file.path(dir_mnp,"CNV_data","CNanalysis4_conumee_REF-F.vh20150715.RData"))

# read clinical data files
DFclinical_full_cohort <- readRDS(pDFclinical_full_cohort)

# read raw Mset: no normalization is performed before CNV analysis 
Mset_raw <- readRDS(pMset_raw)

# select current sentrix from Mset
# sentrixs = list.files(dir_full_cohort,pattern='_Grn.idat')
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

if (Sex == 'F'){
    ref.data = refF.data
}else if (Sex == 'M'){
    ref.data = refM.data
}else{
    stop('stopped in rule CNVplot at sentrix ',sentrix)
}
x <- CNV.fit(cndata, ref.data, anno)
x <- CNV.bin(x)
# x <- CNV.segment(x,alpha = 0.01, nperm = 10000, min.width = 2, undo.splits = "sdundo", undo.SD = 3)


# findMethods('CNV.fit')
# findMethods('CNV.segment')

# Settings Conumee
# x <- CNV.segment(x, alpha = 0.005, nperm = 50000, min.width = 5, undo.splits = "sdundo", undo.SD = 2.2, verbose = 1)
# # Settings DNAcopy
# x <- CNV.segment(y, alpha = 0.01, nperm = 10000, min.width = 2, undo.splits = "sdundo", undo.SD = 3, verbose = 1)
# settings QDNAseq
# if (bin==100) {SDundo=0.10; alph=1e-20} # default = 1e-20 # 0.01 used for PELLL_FS8_a0.01
# QCN.fcnsds <- segmentBins(QCN.fcnsd[,QCN.fcnsd$used.reads > min_used_reads ], undo.splits='sdundo', undo.SD=SDundo, alpha=alph, transformFun="sqrt") # gives 'Performing segmentation: NA
# thus
# x <- CNV.segment(x, alpha = 1e-20, nperm = 10000, min.width = 2, undo.splits = "sdundo", undo.SD = 0.10, verbose = 1, transformFun='sqrt') # variance square root transformation not included
x <- CNV.segment(x, alpha = 1e-20, nperm = 10000, min.width = 2, undo.splits = "sdundo", undo.SD = 0.10, verbose = 1)
cnv <- x

#plot(cnv@bin$ratio-cnv@bin$shift)
# findMethods('CNV.bin')
message('saving cnv in ', pcnv)
saveRDS(cnv, file = pcnv)


# sentrix <- sentrix
# sMset <- Mset_raw[,sentrix]
# cnsentrix <- CNV.load(sMset)
# mappedsMset <- mapToGenome(sMset)
# anno <- CNV.create_anno(array_type="EPIC")
# anno@probes <- subsetByOverlaps(anno@probes, granges(mappedsMset))
# p <- names(anno@probes)
# sMset <- sMset[p,1]
# cnsentrix <- CNV.load(sMset)
# # inside CNV.fit: https://rdrr.io/bioc/conumee/src/R/process.R
# # r <- cor(query@intensity[p, ], ref@intensity[p, ])[1, ] < 0.99
# # r <- cor(cnsentrix@intensity[p, ], cnref@intensity[p, ])[1, ] < 0.99 cor matrix is 1x1, so need more refs samples?
# 
# ref <- '205061430174_R03C01' # 13752_4:LSS24
# refMset <- Mset_raw[p,ref]
# cnref <- CNV.load(refMset)
# # refMset <- mapToGenome(refMset)
# 
# x <- CNV.fit(cnsentrix, cnref, anno)

Sentrix_ID <- sentrix
# FullSample <- DFclinical_full_cohort$FullSample[DFclinical_full_cohort$Sentrix_ID == sentrix]
LSS <- DFclinical_full_cohort$ID[DFclinical_full_cohort$Sentrix_ID == sentrix]
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


# make QDNAseq usable
# q  = readRDS('./../../output-CGH/100kbp/data/100kbp-segmented.rds')
# names(attributes(q))
# attributes(q)
# pData(attributes(q)$phenoData)
# attributes(q)$annotation
# attributes(q)


sink()
