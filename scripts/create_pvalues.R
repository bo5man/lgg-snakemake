###########
# Erik Bosch 
###########
msg <- snakemake@params[["suppressMessages"]]
if (msg){
    suppressMessages(library(minfi))
    suppressMessages(library(openxlsx))
} else{
    library(minfi)
    library(openxlsx)
}

##################
# input file paths
pRGset_mnp_filtered <- snakemake@input[["RGset_mnp_filtered"]]
ptargets_df <- snakemake@input[["targets_df"]]

# output file paths
ppvalues <- snakemake@output[["pvalues"]]
ptargets_df <- snakemake@output[["targets_df"]]

# parameters
test_pvalues <- snakemake@params[['test_pvalues']]
# if(!dir.exists(dir_pca))(dir.create(dir_pca))

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

# read input
RGset_mnp_filtered <- readRDS(pRGset_mnp_filtered)
targets_df <- readRDS(ptargets_df)

pvalues <- detectionP(RGset, type = "m+u")
message('saving pvalues in ppvalues')
saveRDS(pvalues, file = ppvalues)

barplot(colMeans(pvalues),log='y', las=2, cex.names=0.7, ylab="Mean detection p-values")
abline(h=c(0.01,0.05),col="red")
text(1,0.011,'p=0.01',col="red")
text(1,0.055,'p=0.05',col="red")

keepPvalues1 <- colMeans(pvalues) < 0.01 #test_pvalues[0]
keepPvalues5 <- colMeans(pvalues) < 0.05 #test_pvalues[1]
df$SufficientPValues_1 <- keepPvalues1
df$SufficientPValues_5 <- keepPvalues5

df <- df[,c("Sentrix_ID", "FullSample", "SufficientDNA", "SufficientPValues_1", "SufficientPValues_5" "Survival", "Sex", "Age", "AgeMedian", "PFSmonths", "PFSinterval", "Origin")]

message('saving df of targets for clinical data in targets_df')
saveRDS(df, file = ptargets_df)

sink()
sink()
