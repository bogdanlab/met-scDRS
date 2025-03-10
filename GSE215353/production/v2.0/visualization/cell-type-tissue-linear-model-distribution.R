### cell-type-tissue-linear-model-distribution.R ##################################################
# purpose: plot out the distribution computed from brain-region-linear-model.R

### PREAMBLE ######################################################################################
# load in libraries:
# gather the different regions
require(docopt)
require(circlize);
require(ComplexHeatmap);
require(ggplot2);
require(tidyverse);
require(ggrepel);

### VISUALIZE #####################################################################################
# define parameters:
output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
top.visualize <- 10;
system.date <- Sys.Date();

# read in data
args <- commandArgs(trailingOnly = TRUE)

# Extract the name and age
pval.file <- args[1]
long.p.val <- read.csv(file = pval.file)

# subset to the publication traits:
publication.traits <- c(
    'UKB_460K.blood_RBC_DISTRIB_WIDTH',
    'UKB_460K.blood_MONOCYTE_COUNT',
    'UKB_460K.blood_LYMPHOCYTE_COUNT',
    'PASS_Rheumatoid_Arthritis',
    'PASS_Multiple_sclerosis',
    'PASS_IBD_deLange2017',
    'UKB_460K.disease_ASTHMA_DIAGNOSED',
    'UKB_460K.disease_HYPOTHYROIDISM_SELF_REP',
    'UKB_460K.disease_AID_ALL',
    'PASS_Schizophrenia_Pardinas2018',
    'PASS_MDD_Howard2019',
    'PASS_BIP_Mullins2021',
    'UKB_460K.cov_EDU_COLLEGE',
    'UKB_460K.body_BMIz',
    'UKB_460K.cov_SMOKING_STATUS',
    'UKB_460K.biochemistry_Triglycerides',
    'UKB_460K.biochemistry_Testosterone_Male',
    'UKB_460K.body_HEIGHTz',
    'UKB_460K.bmd_HEEL_TSCOREz',
    'UKB_460K.bp_SYSTOLICadjMEDz',
    'PASS_Type_2_Diabetes',
    'UKB_460K.biochemistry_Glucose'
    );
tissue <- c(
    rep('Blood/immune', 9),
    rep('Brain', 6),
    rep('Others', 7)
    );
names(tissue) <- publication.traits;

# subset to the set of publication traits:
plot.mask <- long.p.val$disease %in% publication.traits
long.p.val <- long.p.val[plot.mask, ]
long.p.val$disease_type <- tissue[match(long.p.val$disease, names(tissue))]

# create plot.df:
long.p.val$negative_log10_p <- -log10(long.p.val$bf_p + 1e-323)
long.p.val$annotation <- paste0(long.p.val$cell_type, ':', long.p.val$tissue, ':', long.p.val$disease);

# create a mask for annotation:
mask <- sort(long.p.val$negative_log10_p, decreasing = TRUE)[top.visualize]
long.p.val$annotation <- gsub('Brodmann..1909..', '', long.p.val$annotation)
long.p.val$annotation <- gsub('UKB_460K.', '', long.p.val$annotation)
long.p.val$annotation <- gsub('PASS_', '', long.p.val$annotation)

# Create plot:
gplot <- ggplot(
    long.p.val,
    aes(y = negative_log10_p, x = effect, color = disease_type)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(
        data = subset(long.p.val, negative_log10_p >= mask),
        aes(
            label = annotation,
            size = 1
            ),
        point.padding = 3,
        max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
        nudge_x = 2,
        nudge_y = 2,
        )
    theme_classic() +
    labs(color = "trait") +
    ylab('- log10 Bonforroni correct p value') +
    xlab('effect size') +
    theme(text = element_text(size = 20))

png(
    filename = paste0(output.path, system.date, '-effect-size-tissue-celltype-disease-volcanic-plot-labelled.png'),
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot);
dev.off();

# Create plot:
gplot <- ggplot(
    long.p.val,
    aes(y = negative_log10_p, x = effect, color = disease_type)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    labs(color = "trait") +
    ylab('- log10 Bonforroni correct p value') +
    xlab('effect size') +
    theme(text = element_text(size = 20))

png(
    filename = paste0(output.path, system.date, '-effect-size-tissue-celltype-disease-volcanic-plot.png'),
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot);
dev.off();
