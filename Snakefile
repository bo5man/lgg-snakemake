# Set config dependent parameters
import os

configfile: "config.yaml"
# report: "report/workflow.rst"

DIR_IDATS   = os.path.join(config['paths']['all_idats'], '')
#DIR_IDATS   = os.path.join(config['paths']['selected_idats'], '')
DIR_FULL_COHORT   = os.path.join(config['paths']['full_cohort'], '')

DIR_OUT     = os.path.join(config['paths']['out'], '')
DIR_LOG     = os.path.join(config['paths']['log'], '')
DIR_MNP     = os.path.join(config['paths']['mnp'], '')
DIR_CNV     = os.path.join(DIR_OUT, 'cnv', '')
DIR_CNV_RDS = os.path.join(DIR_CNV, 'RDS', '')

DIR_BETAS   = os.path.join(DIR_OUT, 'betas', '')
DIR_MSET    = os.path.join(DIR_OUT, 'Mset', '')
DIR_MSET_SINGLES    = os.path.join(DIR_OUT, 'Mset', 'singles', '')
DIR_RGSET   = os.path.join(DIR_OUT, 'RGset', '')
DIR_RGSET_SINGLES   = os.path.join(DIR_OUT, 'RGset', 'singles', '')
DIR_QC      = os.path.join(DIR_OUT, 'QC', '')

DIR_TSNE    = os.path.join(DIR_OUT, 'tSNE', '')
DIR_PCA     = os.path.join(DIR_TSNE, 'pca', '')
DIR_RDS     = os.path.join(DIR_TSNE, 'RDS', '')
#DIR_PERPLEXITY  = os.path.join(DIR_TSNE + '005', '')
#PERPLEXITY  = 5  # '5'.zfill(3) = '005' pad 5 with zeros untill 3 characters

DIR_ANALYSIS_PROBES = os.path.join(DIR_OUT, 'analysis_probes', '')
DIR_ANALYSIS_PROBES_RF = os.path.join(DIR_OUT, 'analysis_probes_rf', '')
DIR_ANALYSIS_PROBES_LR = os.path.join(DIR_OUT, 'analysis_probes_lr', '')
DIR_gCIMP = os.path.join(DIR_ANALYSIS_PROBES,'gCIMP_7probes','')
DIR_GLASS = os.path.join(DIR_ANALYSIS_PROBES,'GLASS','')

#PVALUES     = os.path.join(DIR_OUT, 'pvalues.RDS')

sentrixs, = glob_wildcards(DIR_IDATS + '{sentrix}_Grn.idat')
sentrixs_cohort, = glob_wildcards(DIR_FULL_COHORT + '{sentrix_cohort}_Grn.idat')

def getsentrixs():
    SENTRIXS=dict()
    for sentrix in sentrixs:
        sentrixfile=DIR_IDATS+sentrix+"_Grn.idat"
        SENTRIXS[sentrix]=sentrixfile
    return(SENTRIXS)
SENTRIXS = getsentrixs()
#print(SENTRIXS)
def getfullcohortsentrixs():
    SENTRIXS=dict()
    for sentrix in sentrixs_cohort:
        sentrixfile=DIR_FULL_COHORT+sentrix+"_Grn.idat"
        SENTRIXS[sentrix]=sentrixfile
    return(SENTRIXS)
SENTRIXS_COHORT = getfullcohortsentrixs()
#print(SENTRIXS_COHORT)

rule all: #uncomment which branch you want for testing
    input:
        ###### create untill create_singles #####
        #RGsets =     expand(DIR_RGSET_SINGLES + '{sentrix}_RGset.RDS', sentrix = SENTRIXS.keys()),                # create_singles
        #Msets =      expand(DIR_MSET_SINGLES + '{sentrix}_Mset.RDS', sentrix = SENTRIXS.keys()),                  # create_singles
        #Msets_mnp =  expand(DIR_MSET_SINGLES + '{sentrix}_Mset_mnp_filtered.RDS', sentrix = SENTRIXS.keys()),     # create_singles
        #Msets_raw =  expand(DIR_MSET_SINGLES + '{sentrix}_Mset_raw.RDS', sentrix = SENTRIXS.keys()),              # create_singles
        #Msets_noob = expand(DIR_MSET_SINGLES + '{sentrix}_Mset_noob.RDS', sentrix = SENTRIXS.keys()),             # create_singles
        ###### create until rule create_betas  #####
        RGset  = DIR_RGSET + 'RGset.RDS',                             # create_RGset
        Mset  = DIR_MSET + 'Mset.RDS',                             # create_Msets
        Mset_mnp_filtered  = DIR_MSET + 'Mset_mnp_filtered.RDS',   # create_Msets
        Mset_raw    = DIR_MSET + 'Mset_raw.RDS',                   # create_Msets
        Mset_noob    = DIR_MSET + 'Mset_noob.RDS',                   # create_Msets
        betas   = DIR_BETAS + 'betas.RDS',                            # create_betas
        betas_mnp_filtered  = DIR_BETAS + 'betas_mnp_filtered.RDS',   # create_betas
        ##### create untill clinicalDF  #####
        DFclinical_full_cohort =    DIR_OUT + 'DFclinical_full_cohort.RDS',     # create_clinicalDF
        DFclinical_full_gliomas =   DIR_OUT + 'DFclinical_full_gliomas.RDS',    # create_clinicalDF
        DFclinical_full_inhouse =   DIR_OUT + 'DFclinical_full_inhouse.RDS',    # create_clinicalDF
        DFclinical_full_cohort_csv =    DIR_OUT + 'DFclinical_full_cohort.csv',     # create_clinicalDF
        DFclinical_full_inhouse_csv =    DIR_OUT + 'DFclinical_full_inhouse.csv',     # create_clinicalDF
        ##### create untill check_QC  #####
        qcMset_full_cohort =            DIR_QC + 'qcMset_full_cohort.RDS',                  # create_QC
        qcplotMset_full_cohort =        DIR_QC + 'qcplotMset_full_cohort.png',              # create_QC
        densityplotMset_full_cohort =   DIR_QC + 'densityplotMset_full_cohort.png',         # create_QC
        DFclinical_cohort =     DIR_OUT + 'DFclinical_cohort.RDS',  # check_QC
        DFclinical_gliomas =    DIR_OUT + 'DFclinical_gliomas.RDS', # check_QC
        DFclinical_inhouse =    DIR_OUT + 'DFclinical_inhouse.RDS', # check_QC
        DFclinical_cohort_csv = DIR_OUT + 'DFclinical_cohort.csv',     # check_QC
        DFclinical_inhouse_csv =    DIR_OUT + 'DFclinical_inhouse.csv',     # check_QC
        betas_cohort =      DIR_BETAS + 'betas_cohort.RDS',             # check_QC 
        betas_gliomas =     DIR_BETAS + 'betas_gliomas.RDS',            # check_QC 
        betas_inhouse =     DIR_BETAS + 'betas_inhouse.RDS',            # check_QC 
        ##### create untill tSNE  #####
        #pca_cohort =   DIR_PCA + 'pca_cohort.RDS',            # create_tSNE
        #pca_gliomas =       DIR_PCA + 'pca_gliomas.RDS',    # create_tSNE
        pca_inhouse =       DIR_PCA + 'pca_inhouse.RDS',    # create_tSNE
        #tsne_cohort =          DIR_RDS + 'tsne_cohort.RDS',   # create_tSNE
        #tsne_gliomas =              DIR_RDS + 'tsne_gliomas.RDS',       # create_tSNE
        tsne_inhouse =              DIR_RDS + 'tsne_inhouse.RDS',       # create_tSNE
        #tsneplot_cohort_Type =             DIR_TSNE + 'tsneplot_cohort_Type.png',         # create_tSNE
        #tsneplot_cohort_TypeSurvival =     DIR_TSNE + 'tsneplot_cohort_TypeSurvival.png', # create_tSNE
        #tsneplot_cohort_ID =              DIR_TSNE + 'tsneplot_cohort_ID.png',          # create_tSNE
        #tsneplot_gliomas_Type =                 DIR_TSNE + 'tsneplot_gliomas_Type.png',             # create_tSNE
        #tsneplot_gliomas_TypeSurvival =         DIR_TSNE + 'tsneplot_gliomas_TypeSurvival.png',     # create_tSNE
        #tsneplot_gliomas_ID =                  DIR_TSNE + 'tsneplot_gliomas_ID.png',              # create_tSNE
        tsneplot_inhouse_Type =                 DIR_TSNE + 'tsneplot_inhouse_Type.png',             # create_tSNE
        tsneplot_inhouse_TypeSurvival =         DIR_TSNE + 'tsneplot_inhouse_TypeSurvival.png',     # create_tSNE
        tsneplot_inhouse_ID =                  DIR_TSNE + 'tsneplot_inhouse_ID.png',              # create_tSNE
        tsneplot_inhouse_BadSample =                  DIR_TSNE + 'tsneplot_inhouse_BadSample.png',              # create_tSNE
        ####### CNV branch #####
        cnv = expand(DIR_CNV_RDS + '{sentrix_cohort}_cnv.RDS', sentrix_cohort = SENTRIXS_COHORT.keys()),                            # create_CNVplot
        cnvplot = expand(DIR_CNV + '{sentrix_cohort}_cnv.png', sentrix_cohort = SENTRIXS_COHORT.keys()),                            # create_CNVplot
        ###### analysis_probes branch #####
        betas_gCIMP =               DIR_gCIMP + 'betas_gCIMP.RDS',   # analysis_probes
        betas_gCIMP_short =         DIR_gCIMP + 'betas_gCIMP_short.RDS',   # analysis_probes
        betas_gCIMP_long =          DIR_gCIMP + 'betas_gCIMP_long.RDS',   # analysis_probes
        #betas_gCIMP_thresholded  =  DIR_gCIMP + 'betas_gCIMP_thresholded.RDS',   # analysis_probes
        ###TO-DO:betas_gCIMP_short_thresholded = dir_gCIMP + 'betas_gCIMP_short_thresholded.RDS',   # analysis_probes
        ###TO-DO:betas_gCIMP_long_thresholded = dir_gCIMP + 'betas_gCIMP_long_thresholded.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes =       DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes_short = DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes_short.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes_long =  DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes_long.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes =          DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes_short.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes_long.RDS',   # analysis_probes
        heatmap_gCIMP_cohort =              DIR_gCIMP + 'heatmap_gCIMP_cohort.png', # analysis_probes
        heatmap_gCIMP_cohort_short =        DIR_gCIMP + 'heatmap_gCIMP_cohort_short.png', # analysis_probes
        heatmap_gCIMP_cohort_long =         DIR_gCIMP + 'heatmap_gCIMP_cohort_long.png', # analysis_probes
        #heatmap_gCIMP_cohort_thresholded =  DIR_gCIMP + 'heatmap_gCIMP_cohort_thresholded.png', # analysis_probes
        heatmap_treatment_related_cohort =          DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort.png', # analysis_probes
        heatmap_treatment_related_cohort_short =    DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_short.png', # analysis_probes
        heatmap_treatment_related_cohort_long =     DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_long.png', # analysis_probes
        heatmap_hypermodulator_cohort =          DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort.png', # analysis_probes
        heatmap_hypermodulator_cohort_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_short.png', # analysis_probes
        heatmap_hypermodulator_cohort_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_long.png', # analysis_probes
        heatmapRDS_gCIMP_cohort =              DIR_gCIMP + 'heatmap_gCIMP_cohort.RDS', # analysis_probes
        heatmapRDS_gCIMP_cohort_short =        DIR_gCIMP + 'heatmap_gCIMP_cohort_short.RDS', # analysis_probes
        heatmapRDS_gCIMP_cohort_long =         DIR_gCIMP + 'heatmap_gCIMP_cohort_long.RDS', # analysis_probes
        #heatmapRDS_gCIMP_cohort_thresholded =  DIR_gCIMP + 'heatmap_gCIMP_cohort_thresholded.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort =          DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort_short =    DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_short.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort_long =     DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_long.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort =          DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_short.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_long.RDS', # analysis_probes
        ####### analysis_probes_rf branch
        ###TO-DO: DIR_ANALYSIS_PROBES + 'DFfeatures_full_cohort.RDS',  # analysis_probes_rf
        ###TO-DO: DIR_ANALYSIS_PROBES + 'forest_full_cohort.RDS',     # analysis_probes_rf
        ###TO-DO: DIR_ANALYSIS_PROBES + 'forest_roc_full_cohort.RDS',     # analysis_probes_rf
        ####### analysis_probes_lr branch
        glm_fullstats_gCIMP =   DIR_ANALYSIS_PROBES_LR + 'glm_fullstats_gCIMP.RDS', # analysis_probes_lr
        glm_stats_gCIMP =   DIR_ANALYSIS_PROBES_LR + 'glm_stats_gCIMP.xlsx', # analysis_probes_lr
        DFfeatures_gCIMP =  DIR_ANALYSIS_PROBES_LR + 'DFfeatures_gCIMP.RDS',    # analysis_probes_lr
        glm_fullstats_glass_treatment_related =   DIR_ANALYSIS_PROBES_LR + 'glm_fullstats_glass_treatment_related.RDS', # analysis_probes_lr
        glm_stats_glass_treatment_related =   DIR_ANALYSIS_PROBES_LR + 'glm_stats_glass_treatment_related.xlsx', # analysis_probes_lr
        DFfeatures_glass_treatment_related =  DIR_ANALYSIS_PROBES_LR + 'DFfeatures_glass_treatment_related.RDS',    # analysis_probes_lr

rule create_RGset:
    input:
        # sentrix  = DIR_IDATS + '{sentrix}_Grn.idat',
        sentrix  = expand(DIR_IDATS + '{sentrix}_Grn.idat', sentrix = SENTRIXS.keys()),
    output:
        RGset  = DIR_RGSET + 'RGset.RDS',                             # create_RGset
        # RGset  = report(DIR_RGSET + 'RGset.RDS', category="Create raw data sets", subcategory="RedGreen Methylation Channels"),                             # create_RGset
    params:
        dir_rgset   = DIR_RGSET,
        dir_mnp = DIR_MNP,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_RGset/log.log"
    script:
        'scripts/create_RGset.R'

rule create_Msets:
    input:
        RGset  = DIR_RGSET + 'RGset.RDS',                             # create_RGset
    output:
        # Mset  = report(DIR_MSET + 'Mset.RDS', category="Create raw data sets", subcategory="Methylation set"),                             # create_Msets
        Mset  = DIR_MSET + 'Mset.RDS',                             # create_Msets
        Mset_mnp_filtered  = DIR_MSET + 'Mset_mnp_filtered.RDS',   # create_Msets
        Mset_raw    = DIR_MSET + 'Mset_raw.RDS',                   # create_Msets
        # Mset_raw    = report(DIR_MSET + 'Mset_raw.RDS', category="Create raw data sets", subcategory="raw preprocessed Methylation set"),                   # create_Msets
        Mset_noob    = DIR_MSET + 'Mset_noob.RDS',                   # create_Msets
    params:
        dir_mnp = DIR_MNP,
        dir_mset    = DIR_MSET,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_Msets/log.log"
    script:
        'scripts/create_Msets.R'

rule create_betas:
    input:
        Mset  = DIR_MSET + 'Mset.RDS',                             # create_Msets
        Mset_mnp_filtered  = DIR_MSET + 'Mset_mnp_filtered.RDS',   # create_Msets
    output:
        # betas   = report(DIR_BETAS + 'betas.RDS', category="Create raw data sets", subcategory="beta values of Methylation set"),                            # create_betas
        betas   = DIR_BETAS + 'betas.RDS',                            # create_betas
        betas_mnp_filtered  = DIR_BETAS + 'betas_mnp_filtered.RDS',   # create_betas
        # betas_mnp_filtered_sd  = DIR_BETAS + 'betas_mnp_filtered_sd.RDS',     # create_betas_sorted_sd
    params:
        # getBeta_offset  = 100
        # sd_filter = 32000,
        dir_mnp = DIR_MNP,
        dir_betas   = DIR_BETAS,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_betas/log.log"
    script:
        'scripts/create_betas.R'

rule create_singles:
    input:
        sentrix  = DIR_IDATS + '{sentrix}_Grn.idat',
        RGset  = DIR_RGSET + 'RGset.RDS',                             # create_RGset
        Mset  = DIR_MSET + 'Mset.RDS',                             # create_Msets
        Mset_mnp_filtered  = DIR_MSET + 'Mset_mnp_filtered.RDS',   # create_Msets
        Mset_raw    = DIR_MSET + 'Mset_raw.RDS',                   # create_Msets
        Mset_noob    = DIR_MSET + 'Mset_noob.RDS',                   # create_Msets
    output:
        singleRGset  = DIR_RGSET_SINGLES + '{sentrix}_RGset.RDS',                             # create_singles
        singleMset  = DIR_MSET_SINGLES + '{sentrix}_Mset.RDS',                             # create_singles
        singleMset_mnp_filtered  = DIR_MSET_SINGLES + '{sentrix}_Mset_mnp_filtered.RDS',   # create_singles
        singleMset_raw    = DIR_MSET_SINGLES + '{sentrix}_Mset_raw.RDS',                   # create_singles
        singleMset_noob    = DIR_MSET_SINGLES + '{sentrix}_Mset_noob.RDS',                   # create_singles
    params:
        dir_idats   = DIR_IDATS,
        dir_rgset     = DIR_RGSET,
        dir_rgset_singles     = DIR_RGSET_SINGLES,
        dir_mset     = DIR_MSET,
        dir_mset_singles     = DIR_MSET_SINGLES,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_singles/{sentrix}.log"
    script:
        'scripts/create_singles.R'

rule create_clinicalDF:
    input:
        #sentrix  = expand(DIR_FULL_COHORT + '{sentrix_cohort}_Grn.idat', sentrix_cohort=SENTRIXS_COHORT.keys()),
        Mset   = DIR_MSET + 'Mset_mnp_filtered.RDS',                            # create_Mset
        betas   = DIR_BETAS + 'betas_mnp_filtered.RDS',                            # create_betas
    output:
        DFclinical_full_cohort =    DIR_OUT + 'DFclinical_full_cohort.RDS',     # create_clinicalDF
        DFclinical_full_gliomas =   DIR_OUT + 'DFclinical_full_gliomas.RDS',    # create_clinicalDF
        DFclinical_full_inhouse =   DIR_OUT + 'DFclinical_full_inhouse.RDS',    # create_clinicalDF
        # DFclinical_full_inhouse =   report(DIR_OUT + 'DFclinical_full_inhouse.RDS', category="Clinical data sets", subcategory="full data set")    # create_clinicalDF
        DFclinical_full_cohort_csv =    report(DIR_OUT + 'DFclinical_full_cohort.csv', category="Clinical data sets", subcategory="full data set"),     # create_clinicalDF
        DFclinical_full_inhouse_csv =    report(DIR_OUT + 'DFclinical_full_inhouse.csv', category="Clinical data sets", subcategory="full inhouse data set"),     # create_clinicalDF
        Mset_full   = DIR_MSET + 'Mset_full.RDS',                            # create_clinicalDF
        betas_full   = DIR_BETAS + 'betas_full.RDS',                            # create_clinicalDF
    params:
        dir_betas = DIR_BETAS,
        dir_mnp = DIR_MNP,
        inhouse = config['paths']['inhouse'],
        overview = config['paths']['overview'],
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_clinicalDF/log.log"
    script:
        'scripts/create_clinicalDF.R'

rule check_QC:
    input:
        # sentrix  = expand(DIR_FULL_COHORT + '{sentrix_cohort}_Grn.idat', sentrix_cohort=SENTRIXS_COHORT.keys()),
        # RGset  = DIR_RGSET + 'RGset.RDS',                             # create_RGset
        Mset_full   = DIR_MSET + 'Mset_full.RDS',                            # create_clinicalDF
        betas_full   = DIR_BETAS + 'betas_full.RDS',                            # create_clinicalDF
        DFclinical_full_cohort =    DIR_OUT + 'DFclinical_full_cohort.RDS',     # create_clinicalDF
        DFclinical_full_gliomas =   DIR_OUT + 'DFclinical_full_gliomas.RDS',    # create_clinicalDF
        DFclinical_full_inhouse =   DIR_OUT + 'DFclinical_full_inhouse.RDS',    # create_clinicalDF
    output:
        qcMset_full_cohort =            DIR_QC + 'qcMset_full_cohort.RDS',                  # create_QC
        qcplotMset_full_cohort =        report(DIR_QC + 'qcplotMset_full_cohort.png', category="Quality Control"),              # create_QC
        densityplotMset_full_cohort =   report(DIR_QC + 'densityplotMset_full_cohort.png', category="Quality Control"),         # create_QC
        DFclinical_cohort =             DIR_OUT + 'DFclinical_cohort.RDS',  # check_QC
        DFclinical_gliomas =            DIR_OUT + 'DFclinical_gliomas.RDS', # check_QC
        DFclinical_inhouse =            DIR_OUT + 'DFclinical_inhouse.RDS', # check_QC
        DFclinical_cohort_csv =    report(DIR_OUT + 'DFclinical_cohort.csv', category="Clinical data sets", subcategory="Quality Controled data set"),     # check_QC
        DFclinical_inhouse_csv =    report(DIR_OUT + 'DFclinical_inhouse.csv', category="Clinical data sets", subcategory="Quality Controled data set of inhouse sample"),     # check_QC
        # DFclinical_inhouse =            report(DIR_OUT + 'DFclinical_inhouse.RDS', category="Clinical data sets", subcategory="Quality Controled data set"), # check_QC
        betas_cohort =                  DIR_BETAS + 'betas_cohort.RDS',             # check_QC 
        betas_gliomas =                 DIR_BETAS + 'betas_gliomas.RDS',            # check_QC 
        betas_inhouse =                 DIR_BETAS + 'betas_inhouse.RDS',            # check_QC 
        # betas_inhouse =                 report(DIR_BETAS + 'betas_inhouse.RDS', category="Create raw data sets", subcategory="cohort and inhouse sample beta values"),            # check_QC 
    params:
        overview    = config['paths']['overview'],
        dir_full_cohort   = DIR_FULL_COHORT,
        dir_qc     = DIR_QC,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "check_QC/log.log"
    script:
        'scripts/check_QC.R'

rule create_tSNE:
    input:
        betas_full =                 DIR_BETAS + 'betas_full.RDS',
        #betas_cohort =          DIR_BETAS + 'betas_cohort.RDS',             # check_QC 
        #betas_gliomas =         DIR_BETAS + 'betas_gliomas.RDS',            # check_QC 
        betas_inhouse =         DIR_BETAS + 'betas_inhouse.RDS',            # check_QC 
        #DFclinical_cohort =     DIR_OUT + 'DFclinical_cohort.RDS',# create_clinicalDF
        #DFclinical_gliomas =    DIR_OUT + 'DFclinical_gliomas.RDS',    # create_clinicalDF
        DFclinical_inhouse =    DIR_OUT + 'DFclinical_inhouse.RDS',    # create_clinicalDF
        DFclinical_full_inhouse =    DIR_OUT + 'DFclinical_full_inhouse.RDS',    # create_clinicalDF
    output:
        #pca_cohort =        DIR_PCA + 'pca_cohort.RDS',            # create_tSNE
        #pca_gliomas =       DIR_PCA + 'pca_gliomas.RDS',    # create_tSNE
        pca_inhouse =       DIR_PCA + 'pca_inhouse.RDS',    # create_tSNE
        #tsne_cohort =               DIR_RDS + 'tsne_cohort.RDS',        # create_tSNE
        #tsne_gliomas =              DIR_RDS + 'tsne_gliomas.RDS',       # create_tSNE
        #tsne_inhouse =              DIR_RDS + 'tsne_inhouse.RDS',       # create_tSNE
        #tsne_cohort =               report(DIR_RDS + 'tsne_cohort.RDS', category="tSNE", subcategory="cohort"),   # create_tSNE
        #tsne_gliomas =              report(DIR_RDS + 'tsne_gliomas.RDS', category="tSNE", subcategory="cohort and gliomas"),       # create_tSNE
        tsne_inhouse =              DIR_RDS + 'tsne_inhouse.RDS',       # create_tSNE
        #tsneplot_cohort_Type =              report(DIR_TSNE + 'tsneplot_cohort_Type.png', category="tSNE", subcategory="cohort"),         # create_tSNE
        #tsneplot_cohort_TypeSurvival =      report(DIR_TSNE + 'tsneplot_cohort_TypeSurvival.png', category="tSNE", subcategory="cohort"), # create_tSNE
        #tsneplot_cohort_ID =                report(DIR_TSNE + 'tsneplot_cohort_ID.png', category="tSNE", subcategory="cohort"),          # create_tSNE
        #tsneplot_gliomas_Type =             report(DIR_TSNE + 'tsneplot_gliomas_Type.png', category="tSNE", subcategory="cohort and gliomas"),             # create_tSNE
        #tsneplot_gliomas_TypeSurvival =     report(DIR_TSNE + 'tsneplot_gliomas_TypeSurvival.png', category="tSNE", subcategory="cohort and gliomas"),     # create_tSNE
        #tsneplot_gliomas_ID =               report(DIR_TSNE + 'tsneplot_gliomas_ID.png', category="tSNE", subcategory="cohort and gliomas"),              # create_tSNE
        tsneplot_inhouse_Type =             report(DIR_TSNE + 'tsneplot_inhouse_Type.png', category="tSNE", subcategory="cohort and inhouse"),             # create_tSNE
        tsneplot_inhouse_TypeSurvival =     report(DIR_TSNE + 'tsneplot_inhouse_TypeSurvival.png', category="tSNE", subcategory="cohort"),     # create_tSNE
        tsneplot_inhouse_ID =               report(DIR_TSNE + 'tsneplot_inhouse_ID.png', category="tSNE", subcategory="cohort and inhouse"),              # create_tSNE
        tsneplot_inhouse_BadSample =               report(DIR_TSNE + 'tsneplot_inhouse_BadSample.png', category="tSNE", subcategory="cohort and inhouse"),              # create_tSNE
    params:
        # sd_filter = 32000,
        # dir_pca_sd = DIR_PCA_SD,
        dir_betas = DIR_BETAS,
        dir_tsne    = DIR_TSNE,
        dir_pca = DIR_PCA,
        #dir_perplexity = DIR_PERPLEXITY,
        dir_rds = DIR_RDS,
        dir_mset = DIR_MSET,
        dir_mnp = DIR_MNP,
        inhouse = config['paths']['inhouse'],
        overview = config['paths']['overview'],
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_tSNE/log.log"
    script:
        'scripts/create_tSNE.R'

rule create_CNVplot:
    input:
        sentrix  = DIR_FULL_COHORT + '{sentrix_cohort}_Grn.idat',
        Mset_raw    = DIR_MSET + 'Mset_raw.RDS',                   # create_Msets
        DFclinical_full_cohort =    DIR_OUT + 'DFclinical_full_cohort.RDS',# create_clinicalDF
    output:
        cnv =       DIR_CNV_RDS + '{sentrix_cohort}_cnv.RDS',                            # create_CNVplot
        cnvplot =   report(DIR_CNV + '{sentrix_cohort}_cnv.png', category="Results", subcategory="Methylation Copy Number Variation plots"),                            # create_CNVplot
    params:
        dir_full_cohort   = DIR_FULL_COHORT,
        dir_mnp     = DIR_MNP,
        dir_cnv = DIR_CNV,
        dir_cnv_rds = DIR_CNV_RDS,
        # algorithm = ['QDNAseq', 'conumee'],  # make different versions of CNVplots
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "create_CNVplot/{sentrix_cohort}.log"
    script:
        'scripts/create_CNVplot.R'

rule analysis_probes:
    input:
        betas   = DIR_BETAS + 'betas.RDS',                            # create_betas
        DFclinical_cohort =    DIR_OUT + 'DFclinical_cohort.RDS',    # create_tSNE
        DFclinical_gliomas =        DIR_OUT + 'DFclinical_gliomas.RDS',    # create_tSNE
        DFclinical_inhouse =        DIR_OUT + 'DFclinical_inhouse.RDS',    # create_tSNE
    output:
        betas_gCIMP =               DIR_gCIMP + 'betas_gCIMP.RDS',   # analysis_probes
        betas_gCIMP_short =         DIR_gCIMP + 'betas_gCIMP_short.RDS',   # analysis_probes
        betas_gCIMP_long =          DIR_gCIMP + 'betas_gCIMP_long.RDS',   # analysis_probes
        #betas_gCIMP_thresholded  =  DIR_gCIMP + 'betas_gCIMP_thresholded.RDS',   # analysis_probes
        #TO-DO:betas_gCIMP_short_thresholded = dir_gCIMP + 'betas_gCIMP_short_thresholded.RDS',   # analysis_probes
        #TO-DO:betas_gCIMP_long_thresholded = dir_gCIMP + 'betas_gCIMP_long_thresholded.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes =       DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes_short = DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes_short.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes_long =  DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes_long.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes =          DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes_short.RDS',   # analysis_probes
        betas_glass_hypermodulator_342probes_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'betas_hypermodulator_342probes_long.RDS',   # analysis_probes
        heatmap_gCIMP_cohort =              report(DIR_gCIMP + 'heatmap_gCIMP_cohort.png', category="Results", subcategory="Heatmap of selected probes"), # analysis_probes
        heatmap_gCIMP_cohort_short =        DIR_gCIMP + 'heatmap_gCIMP_cohort_short.png', # analysis_probes
        heatmap_gCIMP_cohort_long =         DIR_gCIMP + 'heatmap_gCIMP_cohort_long.png', # analysis_probes
        #heatmap_gCIMP_cohort_thresholded =  DIR_gCIMP + 'heatmap_gCIMP_cohort_thresholded.png', # analysis_probes
        heatmap_treatment_related_cohort =          report(DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort.png', category="Results", subcategory="Heatmap of selected probes"), # analysis_probes
        heatmap_treatment_related_cohort_short =    DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_short.png', # analysis_probes
        heatmap_treatment_related_cohort_long =     DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_long.png', # analysis_probes
        heatmap_hypermodulator_cohort =          report(DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort.png', category="Results", subcategory="Heatmap of selected probes"), # analysis_probes
        heatmap_hypermodulator_cohort_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_short.png', # analysis_probes
        heatmap_hypermodulator_cohort_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_long.png', # analysis_probes
        heatmapRDS_gCIMP_cohort =              DIR_gCIMP + 'heatmap_gCIMP_cohort.RDS', # analysis_probes
        heatmapRDS_gCIMP_cohort_short =        DIR_gCIMP + 'heatmap_gCIMP_cohort_short.RDS', # analysis_probes
        heatmapRDS_gCIMP_cohort_long =         DIR_gCIMP + 'heatmap_gCIMP_cohort_long.RDS', # analysis_probes
        #heatmapRDS_gCIMP_cohort_thresholded =  DIR_gCIMP + 'heatmap_gCIMP_cohort_thresholded.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort =          DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort_short =    DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_short.RDS', # analysis_probes
        heatmapRDS_treatment_related_cohort_long =     DIR_GLASS + 'treatment_related_620probes/' + 'heatmap_treatment_related_cohort_long.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort =          DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort_short =    DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_short.RDS', # analysis_probes
        heatmapRDS_hypermodulator_cohort_long =     DIR_GLASS + 'hypermodulator_342probes/' + 'heatmap_hypermodulator_cohort_long.RDS', # analysis_probes
    params:
        dir_analysis_probes = DIR_ANALYSIS_PROBES,
        dir_gCIMP = DIR_gCIMP,
        dir_glass = DIR_GLASS,
        dir_glass_treatment_related_probes = DIR_GLASS + 'treatment_related_620probes/',
        dir_glass_hypermodulator_probes = DIR_GLASS + 'hypermodulator_342probes/',
        probeset_gCIMP = config['probes']['probeset_gCIMP'],
        probeset_GLASS = config['probes']['probeset_GLASS'],
        probeset_GLASS_treatment_related_sheetname = config['probes']['probeset_GLASS_treatment_related_sheetname'],
        probeset_GLASS_hypermodulator_sheetname = config['probes']['probeset_GLASS_hypermodulator_sheetname'],
        inhouse = config['paths']['inhouse'],
        overview = config['paths']['overview'],
        dir_mnp     = DIR_MNP,
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "analysis_probes/log.log"
    script:
        'scripts/analysis_probes.R'

rule analysis_probes_lr:
    input:
        betas_gCIMP =               DIR_gCIMP + 'betas_gCIMP.RDS',   # analysis_probes
        betas_glass_treatment_related_620probes =       DIR_GLASS + 'treatment_related_620probes/' + 'betas_treatment_related_620probes.RDS',   # analysis_probes
        DFclinical_cohort =    DIR_OUT + 'DFclinical_cohort.RDS',    # create_clinicalDF
    output:
        glm_fullstats_gCIMP =   DIR_ANALYSIS_PROBES_LR + 'glm_fullstats_gCIMP.RDS', # analysis_probes_lr
        glm_stats_gCIMP =   report(DIR_ANALYSIS_PROBES_LR + 'glm_stats_gCIMP.xlsx', category="Results", subcategory="Prediction of long-short survivors"), # analysis_probes_lr
        DFfeatures_gCIMP =  DIR_ANALYSIS_PROBES_LR + 'DFfeatures_gCIMP.RDS',    # analysis_probes_lr
        glm_fullstats_glass_treatment_related =   DIR_ANALYSIS_PROBES_LR + 'glm_fullstats_glass_treatment_related.RDS', # analysis_probes_lr
        glm_stats_glass_treatment_related =   report(DIR_ANALYSIS_PROBES_LR + 'glm_stats_glass_treatment_related.xlsx', category="Results", subcategory="Prediction of long-short survivors"), # analysis_probes_lr
        DFfeatures_glass_treatment_related =  DIR_ANALYSIS_PROBES_LR + 'DFfeatures_glass_treatment_related.RDS',    # analysis_probes_lr
    params:
        DIR_OUT =   DIR_OUT,
        dir_analysis_probes = DIR_ANALYSIS_PROBES,
        dir_analysis_probes_rf = DIR_ANALYSIS_PROBES_RF,
        dir_analysis_probes_lr = DIR_ANALYSIS_PROBES_LR,
        dir_gCIMP = DIR_gCIMP,
        dir_glass = DIR_GLASS,
        dir_glass_treatment_related_probes = DIR_GLASS + 'treatment_related_620probes/',
        dir_glass_hypermodulator_probes = DIR_GLASS + 'hypermodulator_342probes/',
        probeset_gCIMP = config['probes']['probeset_gCIMP'],
        probeset_GLASS = config['probes']['probeset_GLASS'],
        probeset_GLASS_treatment_related_sheetname = config['probes']['probeset_GLASS_treatment_related_sheetname'],
        probeset_GLASS_hypermodulator_sheetname = config['probes']['probeset_GLASS_hypermodulator_sheetname'],
        suppressMessages    = config['options']['suppressMessages'],
    log:
        DIR_LOG + "analysis_probes_lr/log.log"
    script:
        'scripts/analysis_probes_lr.R'

#rule analysis_probes_rf:
## TO-DO: make random forest work
#    input:
#        betas_glass_treatment_related_620probes =       DIR_GLASS + 'betas_treatment_related_620probes.RDS',   # analysis_probes
#        betas_glass_hypermodulator_342probes =          DIR_GLASS + 'betas_hypermodulator_342probes.RDS',   # analysis_probes
#        # DFclinical_gCIMP =                                  DIR_gCIMP + 'DFclinical_gCIMP.RDS',  # analysis_probes
#        DFclinical_glass_treatment_related_620probes =      DIR_GLASS + 'treatment_related_620probes/' + 'DFclinical_treatment_related_620probes.RDS',   # analysis_probes
#        DFclinical_glass_hypermodulator_342probes =         DIR_GLASS + 'hypermodulator_342probes/' + 'DFclinical_hypermodulator_342probes.RDS',   # analysis_probes
#    output:
#        DFfeatures_full_cohort  =                              DIR_ANALYSIS_PROBES + 'DFfeatures_full_cohort.RDS',  # analysis_probes_rf
#        forest_full_cohort =                                   DIR_ANALYSIS_PROBES + 'forest_full_cohort.RDS',     # analysis_probes_rf
#        forest_roc_full_cohort =                               DIR_ANALYSIS_PROBES + 'forest_roc_full_cohort.RDS',     # analysis_probes_rf
#    params:
#        dir_analysis_probes = DIR_ANALYSIS_PROBES,
#        dir_analysis_probes_rf = DIR_ANALYSIS_PROBES_RF,
#        dir_gCIMP = DIR_gCIMP,
#        dir_glass = DIR_GLASS,
#        dir_glass_treatment_related_probes = DIR_GLASS + 'treatment_related_620probes/',
#        dir_glass_hypermodulator_probes = DIR_GLASS + 'hypermodulator_342probes/',
#        probeset_gCIMP = config['probes']['probeset_gCIMP'],
#        probeset_GLASS = config['probes']['probeset_GLASS'],
#        probeset_GLASS_treatment_related_sheetname = config['probes']['probeset_GLASS_treatment_related_sheetname'],
#        probeset_GLASS_hypermodulator_sheetname = config['probes']['probeset_GLASS_hypermodulator_sheetname'],
#        suppressMessages    = config['options']['suppressMessages'],
#    log:
#        DIR_LOG + "analysis_probes_rf/log.log"
#    script:
#        'scripts/analysis_probes_rf.R'
#
#rule CGHtest: # from QDNAseq-snakemake
#    input:
#        RegionsCGH=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-RegionsCGH.rds",
#    output:
#        freqPlotCompare=DIR_OUT + "{ACEbinSize}kbp/CGHtest/freqPlotCompare.png", 
#        CGHtest=DIR_OUT + "{ACEbinSize}kbp/CGHtest/CGHtest.Rdata",
#        plotPFDR=DIR_OUT + "{ACEbinSize}kbp/CGHtest/plotPFDR.png",
#        combined=DIR_OUT + "{ACEbinSize}kbp/CGHtest/combined.png"
#    params:
#        outputdir=DIR_OUT + "{ACEbinSize}kbp/CGHtest/",
#        clinicaldataPath=config["CGHtest"]["clinicaldataPath"],
#        columnSampleNames=config["CGHtest"]["columnSampleNames"],
#        ClassSamples=config["CGHtest"]["ClassSamples"],
#        columnClassSamples=config["CGHtest"]["columnClassSamples"],
#        suppressMessages=config["pipeline"]["suppressMessages"]
#    log:DIR_OUT + DIR_LOG + "{ACEbinSize}kbp/CGHtest_log.tsv"
#    script:
#        "scripts/Run_CGHtest.R"
#
#rule create_pvalues:
#    input:
#        RGset_mnp_filtered  = expand(DIR_RGSET + '{sentrix}_RGset_mnp_filtered.RDS', sentrix = SENTRIXS.keys()),    # create_RGset
#        targets_df  = DIR_OUT + 'targets_df.RDS',                                # create_targets_df
#    output:
#        pvalues = PVALUES,
#        targets_df  = DIR_OUT + 'targets_df.RDS',                                # create_targets_df
#    params:
#        test_pvalues    = [0.01, 0.05]
#        suppressMessages    = config['options']['suppressMessages'],
#    log:
#        DIR_LOG + "create_pvalues.log"
#    script:
#        'scripts/create_pvalues.R'

