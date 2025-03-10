### Visualize MC test.R ###########################################################################
# purpose: aggregate the MC test result from scDRS downstream and then visualize them:

### PREAMBLE ######################################################################################
# specifiy input variables:
require(docopt);
require(circlize);
require(ComplexHeatmap);

'Usage:
    visualize-MC-test.R [--result_dir <result> --field <field> --ordering <order> --out <output>]

Options:
    --result_dir the directories where the result of the group level analysis are located
    --field which of the field would we like to aggregate, can be hetero_mcp or assoc_mcp [default: assoc_mcp]
    --ordering [default: scdrs]
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
group.analysis.dir <- opts$result_dir;
group.index <- opts$field;
output.path <- opts$out;
ordering <- opts$ordering;
system.date <- Sys.Date();

# for testing purposese:
# group.analysis.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v1.2/mch/bimodality-QC-75traits-with-cov/heterogeneity-analysis/';
# group.index <- 'assoc_mcp';
# output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-group-analysis-aggregated', group.index);
# ordering <- 'scdrs'

# load in the data:
group.result.file <- list.files(group.analysis.dir, pattern = 'scdrs_group.annotation');

# initiate place holder:
result.list <- vector('list', length = length(group.result.file));
names(result.list) <- group.result.file;
subset.list <- result.list;

### RESULT EXTRACTION #############################################################################
# extract out result:
for (file in group.result.file) {
    # for each of item in the result list, read in the result:
    result.file <- read.table(
        file = paste0(group.analysis.dir, file),
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    
    result.list[[file]] <- result.file;

    # now extract the data:
    subset.list[[file]] <- result.file[, group.index, drop = FALSE]
    }

# combind the results across traits:
aggregated.result <- Reduce('cbind', subset.list);

# apply fdr correction for each of column:
adjusted.result <- apply(aggregated.result, 2, p.adjust, method = 'fdr');
colnames(adjusted.result) <- gsub('.scdrs_group.annotation', '', group.result.file);

# output the data:
csv.output.path <- paste0(output.path, '.csv');
write.table(
    adjusted.result,
    file = csv.output.path,
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZE #####################################################################################
# visualize the matrix as a heatmap:
# p value correction:
significance.matrix <- t(-log10(adjusted.result));

# check if publication trait is in there:
column.order <- c(
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

if (ordering == 'scdrs'){
    significance.matrix <- significance.matrix[column.order, ]
} else {
    cat('no trait ordering applied, plotting as is\n')
}

# make visualization:
# design olor function:
cutoff <- -log10(0.05);

col.fun <- colorRamp2(
    c(
        0,
        max(significance.matrix, na.rm = TRUE)
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    title_position = 'topcenter',
    direction = 'horizontal',
    at = c(
        0,
        signif(max(significance.matrix, na.rm = TRUE), 3)
        )
    );

# plot out the heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix),
    name = '-log10(P)',
    col = col.fun,
    na_col = 'gray',
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "black", lwd = 2),
    cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(as.matrix(significance.matrix)[i, j]) && as.matrix(significance.matrix)[i, j] > cutoff) {
                grid.text("*", x, y, gp = gpar(col = "black", fontsize = 16))
            }
        },
    width = unit(10 * ncol(significance.matrix),"mm"),
    height = unit(10 * nrow(significance.matrix),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    heatmap_legend_param = heatmap.legend.param
    );

# automatically detect plot size and make plot:
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

heatmap.output.path <- paste0(output.path, '.png')
png(
    filename = heatmap.output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, padding = unit(c(10, 10, 10, 70), "mm"), heatmap_legend_side = "bottom");
dev.off();
