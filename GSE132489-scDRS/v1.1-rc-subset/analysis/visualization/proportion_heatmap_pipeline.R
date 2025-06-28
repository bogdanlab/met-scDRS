### cell-type-proportion-heatmap.R ################################################################
# purpose: plot out the proportion of cells in each cell type that are significant in each traits

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--dir <scdrs> --meta_data <meta> --modality <modality> --p_cutoff <p_cutoff>]

Options:
    --dir directory path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --modality either mch or mcg
    --p_cutoff FDR corrected p value to call significant
]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$dir;
modality <- opts$modality
system.date <- Sys.Date();
p.cutoff <- as.numeric(opts$p_cutoff);

# load libraries:
library(ggplot2);
library(ComplexHeatmap);
library(circlize);
library(readxl);
library(data.table);

# create place holder:
sig.matrix.proportion.collection <- vector('list', length = length(c('mch', 'mcg')));
names(sig.matrix.proportion.collection) <- c('mch', 'mcg');

score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# load meta data:
meta = read.table(file = meta.data.path, sep = ',', header = TRUE, row.names = 1)

# read into the empty list:
for (result in score.files) {
    risk.score[[result]] <- read.table(
        file = result,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# load in a function that would calculate the significant matrix for us:
source('/u/home/l/lixinzhe/project-github/spells/met-scDRS/get-significant-matrix.R')

### ANALYSIS: SIGNIFICANCE CELLS ##################################################################
# take a look at the significance matrix by spatial dissection:
feature.attribute <- c(
    'RegionName',
    'MajorRegion',
    'SubRegion',
    'CellClass',
    'MajorType',
    'SubType'
    );

# get the proportion and count significance matrix:
sig.matrix.proportion <- lapply(
    X = feature.attribute,
    get.significant.matrix,
    risk.score.collection = risk.score,
    meta.data = meta,
    plot.type = 'proportion',
    p.cutoff = p.cutoff
    );
sig.matrix.count <- lapply(
    X = feature.attribute,
    get.significant.matrix,
    risk.score.collection = risk.score,
    meta.data = meta,
    plot.type = 'count',
    p.cutoff = p.cutoff
    );
names(sig.matrix.proportion) <- names(sig.matrix.count) <- feature.attribute;

# define subset of traits that we are interested in:
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

# define space that we will annotate as blood; brain; other traits:
row.split = c(
    rep('Blood/immune', 9),
    rep('Brain', 6),
    rep('Others', 7)
    );

### VISUALIZATION: MajorREGION ####################################################################
# in our data, they are hierarchical as follow: MajorRegion -> SubRegion -> RegionName
annotation <- 'MajorRegion'
significance.matrix <- sig.matrix.proportion[[annotation]];

# plot out the significance heatmap proportion
# we will only keep the features where within each level, there needs to be > 100 cells for plots:
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];
plot.df <- significance.matrix[publication.traits, feature.keep]

# define color function that goes from white to red:
col.fun <- colorRamp2(
    c(0, 0.3, 0.5, 0.7, 1),
    c('white', '#fcae91', '#fb6a4a', '#de2d26', '#de2d26')
    );
heatmap.legend.param <- list(at = c(0, 0.5, 1));

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    # column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

# output the results:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZATION: SubRegion ######################################################################
annotation <- 'SubRegion'
significance.matrix <- sig.matrix.proportion[[annotation]];

# plot out the significance heatmap proportion
# we will only keep the features where within each level, there needs to be > 100 cells for plots:
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];
plot.df <- significance.matrix[publication.traits, feature.keep]

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    # column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

# output the results:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZATION: RegionName #####################################################################
annotation <- 'RegionName'
significance.matrix <- sig.matrix.proportion[[annotation]];

# plot out the significance heatmap proportion
# we will only keep the features where within each level, there needs to be > 100 cells for plots:
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];
plot.df <- significance.matrix[publication.traits, feature.keep]

# define another colour function:
col.fun <- colorRamp2(
    c(0, 1),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(at = c(0, 1));

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    # column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 30,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

# output the results:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZATION: CellClass ##################################################################
# cell identity is hierarchically ordered like so: CellClass -> MajorType -> SubType
annotation <- 'CellClass'
significance.matrix <- sig.matrix.proportion[[annotation]];

# plot out the significance heatmap proportion
# we will only keep the features where within each level, there needs to be > 100 cells for plots:
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];
plot.df <- significance.matrix[publication.traits, feature.keep]

# define color function that goes from white to red:
col.fun <- colorRamp2(
    c(0, 0.3, 0.5, 0.7, 1),
    c('white', '#fcae91', '#fb6a4a', '#de2d26', '#de2d26')
    );
heatmap.legend.param <- list(at = c(0, 0.5, 1));

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    # column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 10,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

# output the results:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZATION - MajorType #####################################################################
annotation <- 'MajorType'
significance.matrix <- sig.matrix.proportion[[annotation]];

# plot out the significance heatmap proportion
# we will only keep the features where within each level, there needs to be > 100 cells for plots:
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];

# define column order for so that it is ordered:
all.region <- c(
    'ODC',
    'OPC',
    'ASC',
    'MGC',
    'EC',
    'PC',
    'VLMC',
    'VLMC-Pia'
    );
ANP.region <- c('ANP');
cortex.region <- c(
    'IT-L23',
    'IT-L4',
    'IT-L5',
    'IT-L6',
    'PT-L5',
    'CT-L6',
    'NP-L6',
    'L6b'
    );
OLF.exc.region <- c('CLA','EP','OLF-Exc');
HIP.region <- c(
    'IG-CA2',
    'CA1',
    'CA3',
    'CA3-St18',
    'Gfra1',
    'DG-po',
    'DG'
    );
cortex.hip.olf.region <- c('MGE-Pvalb', 'MGE-Sst', 'CGE-Vip', 'CGE-Lamp5', 'Unc5c');
CNU.region <- c(
    'Chd7',
    'LSX-Inh',
    'PAL-Inh',
    'Foxp2',
    'MSN-D1',
    'MSN-D2',
    'D1L-PAL'
    );
OLF.region <- c('D1L-Fstl4', 'OLF');
column.order <- c(
    all.region,
    ANP.region,
    cortex.region,
    OLF.exc.region,
    HIP.region,
    cortex.hip.olf.region,
    CNU.region,
    OLF.region
    );

stopifnot(length(column.order) == length(unique(meta[, annotation])));
stopifnot(column.order %in% unique(meta[, annotation]))

column.split <- c(
    rep('ALL', length(all.region)),
    rep('ANP', length(ANP.region)),
    rep('Cortex', length(cortex.region)),
    rep('OLF', length(OLF.exc.region)),
    rep('HIP', length(HIP.region)),
    rep('Cortex;HIP;OLF', length(cortex.hip.olf.region)),
    rep('CNU', length(CNU.region)),
    rep('OLF', length(OLF.region))
    );

plot.df <- significance.matrix[publication.traits, column.order]

# define color function that goes from white to red:
col.fun <- colorRamp2(
    c(0, 0.3, 0.5, 0.7, 1),
    c('white', '#fcae91', '#fb6a4a', '#de2d26', '#de2d26')
    );
heatmap.legend.param <- list(at = c(0, 0.5, 1));

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = column.order,
    width = unit(10 * length(column.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 25,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

# output the results:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZATION - SubType #######################################################################
# for diverse population of cells, we will also want to plot out the sub region:
# number of SubTypes in Major type > 4 would be considered diverse.
annotation <- 'SubType'
significance.matrix <- sig.matrix.proportion[[annotation]];
feature.keep <- names(table(meta[, annotation]))[table(meta[, annotation]) > 100];

# first calculate the number of subtypes within each of the major cell type:
subtypes.in.major <- rowSums(table(meta$MajorType, meta$SubType)[column.order, ] != 0);
diverse.major.types <- names(subtypes.in.major)[subtypes.in.major >= 5];

# get the column order
column.order <- column.split <- {}
for(diverse.major.type in diverse.major.types) {
    # get the cell id within that major type:
    cells.in.diverse.major <- rownames(meta)[meta$MajorType %in% diverse.major.type];
    # get the cell subtype within that major type:
    diverse.minor.types <- unique(meta[cells.in.diverse.major, 'SubType']);

    # concatenate the column order:
    column.order <- c(column.order, diverse.minor.types);
    # obtain the column split:
    column.split <- c(column.split, rep(diverse.major.type, length(diverse.minor.types)));
    }

keep.index <- column.order %in% feature.keep;
column.split <- column.split[keep.index];
column.split <- factor(column.split, levels=unique(column.split))
feature.keep <- column.order[keep.index];

plot.df <- significance.matrix[publication.traits, feature.keep]

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 45,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

### VISUALIZATION FOR PT-L5, ANP, GFRA1, DG-po #################################################### 
annotation <- 'SubType'
significance.matrix <- sig.matrix.proportion[[annotation]];
feature.keep <- c(
    grep('ANP', colnames(significance.matrix)),
    grep('PT-L5', colnames(significance.matrix)),
    grep('Gfra1', colnames(significance.matrix)),
    grep('DG-po', colnames(significance.matrix))
    );
feature.keep <- colnames(significance.matrix)[feature.keep];

column.split <- c(
    rep('ANP', length(grep('ANP', colnames(significance.matrix)))),
    rep('PT-L5', length(grep('PT-L5', colnames(significance.matrix)))),
    rep('Gfra1', length(grep('Gfra1', colnames(significance.matrix)))),
    rep('DG-po', length(grep('DG-po', colnames(significance.matrix))))
    );
column.split <- factor(column.split, levels=unique(column.split))
plot.df <- significance.matrix[publication.traits, feature.keep]

# plot out the significance heatmap:
plot <- Heatmap(
    as.matrix(plot.df),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = feature.keep,
    cluster_columns = FALSE,
    width = unit(10 * length(feature.keep),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0(
    '/u/home/l/lixinzhe/project-geschwind/plot/',
    system.date,
    '-',
    modality,
    '-',
    annotation,
    '-ANP-PTL5-Gfra1-DGpo-gse132489-fraction-74traits-cell-class-significance-proportion',
    '.png'
    );
png(
    filename = output.path,
    width = 45,
    height = 14,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();
