### cell_type_plot.R ##############################################################################
# purpose: plot the cell type umap plot:

### PREAMBLE ######################################################################################
# load in libraries:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--dir <scdrs> --meta_data <meta> --field <group> --out <output> --cutoff <p> --plot_type <count> --color_scale_max <maxima>]

Options:
    --dir directory path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field name in meta which you would like to compute proportion of significant cells in
    --cutoff p value cutoff that user specifies
    --out path to output file
    --plot_type type of plot for output: either count or proportion is accepted [default: count]
    --color_scale_max the maximum color scale for plotting the heatmap. 
]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$dir;
group.index <- opts$field;
output.path <- opts$out;
plot.type <- opts$plot_type;
p.cutoff <- as.numeric(opts$cutoff);
max.color <- as.numeric(opts$color_scale_max)
system.date <- Sys.Date();

library(ggplot2)
library(circlize);
library(Seurat);
library(data.table)
library(ComplexHeatmap)

# define paths:
date <- Sys.Date()

# load in the data:
meta <- read.csv(
    header = TRUE,
    sep = '\t',
    row.names = 1,
    file = meta.data.path
    );

# label the '' as 'unknown'
meta$newL1 = meta$L1
meta$newL2 = meta$L2
meta$newL3 = meta$L3
meta$adjusted_L3 = meta$L3
meta$fine2_age_groups = meta$age_group

# label the '' as 'unknown'
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

# for L2 and L3, add the word Exc and Inh on top:
meta$newL2 = paste0(meta$newL1, '_', meta$newL2)
meta$newL3 = paste0(meta$newL1, '_', meta$newL3)

# read in the data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

mdd = data.frame(
    fread(
        file =  paste0(scDRS.directory, "/PASS_MDD_Howard2019.score.gz"),
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1
    )

# format the meta data:
new_name = gsub('.allc.tsv.gz', '', rownames(mdd))

# get the common cells:
common_cells = intersect(new_name, rownames(meta))
meta = meta[common_cells, ]

# read into the empty list:
for (result in score.files) {
    met_scdrs_result <- data.frame(
        fread(
        file = result,
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1)
    rownames(met_scdrs_result) = new_name
    met_scdrs_result = met_scdrs_result[common_cells, ]
    met_scdrs_result$fdr = p.adjust(met_scdrs_result$pval, method = 'fdr')
    risk.score[[result]] = met_scdrs_result
    }

# simplify names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

### ANALYSIS ######################################################################################
meta[meta[, group.index] == '', group.index] = 'unknown'
cell.type <- unique(meta[, group.index]);
significance.matrix <- matrix(NA, nrow = length(risk.score), ncol = length(cell.type));
rownames(significance.matrix) <- names(risk.score);
colnames(significance.matrix) <- cell.type;

# create the traits by cell type matrix:
for (result in names(risk.score)) {
    # find the index for significant cells:
    significant.cell <- rownames(risk.score[[result]])[
        p.adjust(risk.score[[result]]$pval, method = 'fdr') < p.cutoff
        ];

    # for each cell type, check number of significant cells are part of that cell type:
    for (type in cell.type) {
        # locate cell id that belong in the cell type:
        cell.type.cell <- rownames(meta)[meta[, group.index] %in% type];

        # find the number of cells that are in each of the cell type category:
        if (plot.type == 'count') {
            significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell);
            } else {
                significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell) /
                    length(cell.type.cell);
            }
        }
    }

write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );
    
# grab out some publication traits:
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
    'PASS_Alzheimers_Jansen2019',
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

cell.types = unique(meta[, group.index])
excitatory = sort(cell.types[grep('^Exc', cell.types)])
inhibitory = sort(cell.types[grep('^Inh', cell.types)])
others = setdiff(cell.types, c(excitatory, inhibitory))

cell.type.order = c(
    excitatory,
    inhibitory,
    others
    )

# remove the word pass from the traits:
publication.traits <- gsub('PASS_', '', publication.traits)
rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))

# plot the heatmap using the subset of traits and the ordered cell types:
if ((group.index %in% c('newL1', 'newL2', 'newL3')) & all(publication.traits %in% rownames(significance.matrix))) {

    # first design the color function depending on how many significant cells there are:
    if (plot.type == 'count'){
        if(max(significance.matrix[publication.traits, cell.type.order]) > 100) {
        col.fun <- colorRamp2(
            c(
                0,
                100,
                max(significance.matrix[publication.traits, cell.type.order])
                ),
            c('white', '#fc9272','#de2d26')
            );
        heatmap.legend.param <- list(
            at = c(
                0,
                100,
                round(
                    max(significance.matrix[publication.traits, cell.type.order]),
                    digit = -2 # round to nearest hundreds
                    )
                )
            );
        } else {
        col.fun <- colorRamp2(
            c(
                0,
                max(significance.matrix[publication.traits, cell.type.order])
                ),
            c('white', '#de2d26')
            );
        heatmap.legend.param <- list(
            at = c(
                0,
                max(significance.matrix[publication.traits, cell.type.order])
                )
            );
        }
        } else {
            col.fun <- colorRamp2(c(0, max.color), c('white', '#de2d26'))
            heatmap.legend.param <- list(at = c(0, max.color));
        }

    # plot out the heatmap:
    # define split pattern:
    row.split = c(
        rep('Blood/immune', 9),
        rep('Brain', 7),
        rep('Others', 7)
        );
    column.split = c(
        rep('Excitatory', length(excitatory)),
        rep('Inhibitory', length(inhibitory)),
        rep('Others', length(others))
        )

    plot <- Heatmap(
        as.matrix(significance.matrix)[publication.traits, cell.type.order],
        name = 'Sig. cells',
        col = col.fun,
        rect_gp = gpar(col = "black", lwd = 2),
        row_order = publication.traits,
        column_order = cell.type.order,
        width = unit(10 * length(cell.type.order),"mm"),
        height = unit(10 * length(publication.traits),"mm"),
        column_names_gp = grid::gpar(fontsize = 15),
        row_names_gp = grid::gpar(fontsize = 15),
        row_split = row.split,
        column_split = column.split,
        heatmap_legend_param = heatmap.legend.param
        );
    plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

    # measure the size of the heatmap:
    heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
    heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

    # use the measured width and height for drawing:
    png(
        filename = output.path,
        width = heatmap.width,
        height = heatmap.height,
        units = 'in',
        res = 400
        );
    draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
    dev.off();

    } else {
        publication.traits = intersect(publication.traits, rownames(significance.matrix))
        if (plot.type == 'proportion') {
            # design color function:
            col.fun <- colorRamp2(
                c(
                    0,
                    max.color
                    ),
                c('white', '#de2d26')
                );
            heatmap.legend.param <- list(
                at = c(
                    0,
                    max.color
                    )
                );
            } else {
            # design color function:
            col.fun <- colorRamp2(
                c(
                    0,
                    max(significance.matrix[publication.traits, ])
                    ),
                c('white', '#de2d26')
                );
            heatmap.legend.param <- list(
                at = c(
                    0,
                    max(significance.matrix[publication.traits, ])
                    ),
                c('white', '#de2d26')
                );
            }

        # plot out the heatmap:
        plot <- Heatmap(
            as.matrix(significance.matrix)[publication.traits, cell.type.order],
            name = 'Sig. cells',
            col = col.fun,
            row_order = publication.traits,
            column_order = cell.type.order,
            width = unit(10 * length(cell.type.order),"mm"),
            height = unit(10 * length(publication.traits),"mm"),
            rect_gp = gpar(col = "black", lwd = 2),
            column_names_gp = grid::gpar(fontsize = 15),
            row_names_gp = grid::gpar(fontsize = 15),
            heatmap_legend_param = heatmap.legend.param
            );
        plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

        # measure the size of the heatmap:
        heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
        heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

        # use the measured width and height for drawing:
        png(
            filename = output.path,
            width = heatmap.width,
            height = heatmap.height,
            units = 'in',
            res = 400
            );
        draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
        dev.off();
    }
