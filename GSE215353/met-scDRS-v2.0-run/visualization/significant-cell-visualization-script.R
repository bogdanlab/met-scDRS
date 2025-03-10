### significant-cells-visualization-script.R ######################################################
# purpose: for a directory containing scDRS scores, 
# visualize the number of cells that are significant across all traits
# also visualize the proportion of cell types that are significant given a meta file

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--dir <scdrs> --meta_data <meta> --field <group> --out <output> --cutoff <p> --plot_type <count>]

Options:
    --dir directory path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field name in meta which you would like to compute proportion of significant cells in
    --cutoff p value cutoff that user specifies
    --out path to output file
    --plot_type type of plot for output: either count or proportion is accepted [default: count]

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$dir;
group.index <- opts$field;
output.path <- opts$out;
plot.type <- opts$plot_type;
p.cutoff <- as.numeric(opts$cutoff);
system.date <- Sys.Date();

# function testing:
# meta.data.path <- "/u/project/geschwind/lixinzhe/data/bimodal-asd-visualization-meta.txt"
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/fraction-ASD-methyl/mch/KC/all/"
# group.index <- 'label'
# p.cutoff <- 0.1
# output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
# plot.type <- 'proportion'

# load libraries:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# read data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

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

### ANALYSIS: CELL TYPE NUMBER VISUALIZATION ######################################################
# for each of the result, find out the number of significant cells in each cell type:
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

# output the trait by cell type matrix as supplementary material:
write.table(
    significance.matrix,
    file = gsub('png', 'csv', output.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

# grab out the set of traits that KANGCHENG have used in scDRS publication:
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

cell.type.order <- c(
    "L2/3-IT",
    "L4-IT",
    "L5-ET",
    "L5-IT",
    "L5/6-NP",
    "L6-CT",
    "L6-IT",
    "L6-IT-Car3",
    "L6b",
    "Amy-Exc",
    "CA1",
    "CA3",
    "DG",
    "HIP-Misc1",
    "HIP-Misc2",
    "CB",
    "Chd7",
    "Foxp2",
    "MSN-D1",
    "MSN-D2",
    "PKJ",
    "PN",
    "Lamp5",
    "Lamp5-Lhx6",
    "Pvalb",
    "Pvalb-ChC",
    "Sncg",
    "Sst",
    "SubCtx-Cplx",
    "THM-Exc",
    "THM-Inh",
    "THM-MB",
    "Vip",
    "ASC",
    "EC",
    "MGC",
    "ODC",
    "OPC",
    "PC",
    "VLMC"
    );

# plot the heatmap using the subset of traits and the ordered cell types:
# define color function:
# make the heatmap:
if (all(cell.type.order %in% colnames(significance.matrix)) & all(publication.traits %in% rownames(significance.matrix))) {

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
            col.fun <- colorRamp2(c(0, 1), c('white', '#de2d26'))
            heatmap.legend.param <- list(at = c(0, 1));
        }

    # plot out the heatmap:
    # define split pattern:
    row.split = c(
        rep('Blood/immune', 9),
        rep('Brain', 6),
        rep('Others', 7)
        );
    column.split = c(
        rep('Excitatory', 15),
        rep('Inhibitory', 18),
        rep('Others', 7)
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
        if (plot.type == 'proportion') {
            # design color function:
            col.fun <- colorRamp2(
                c(
                    0,
                    1
                    ),
                c('white', '#de2d26')
                );
            heatmap.legend.param <- list(
                at = c(
                    0,
                    1
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
            as.matrix(significance.matrix)[publication.traits, ],
            name = 'Sig. cells',
            col = col.fun,
            row_order = publication.traits,
            rect_gp = gpar(col = "black", lwd = 2),
            width = unit(5 * ncol(significance.matrix),"mm"),
            height = unit(5 * nrow(significance.matrix),"mm"),
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
