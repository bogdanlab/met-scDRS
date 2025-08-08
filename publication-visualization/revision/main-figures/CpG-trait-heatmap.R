### brain-trait-heatmap.R #########################################################################
# purpose: draw a heatmap that looks at the proportion of significant cells per cell type

### PREAMBLE ######################################################################################
# load packages:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# load in data:
significance.matrix <- read.table(
    file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-08-07-v1.1rc-GSE215353-subset-fraction-mcg-cell-type-significance-proportion.csv',
    sep = ',',
    row.names = 1,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
    );

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

system.date <- Sys.Date()

### VISUALIZATION #################################################################################
# find out the set of brain traits:
trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# select traits to plot:
publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# reformat traits:
publication.traits <- gsub('PASS_', '', publication.traits)
publication.traits <- gsub('UKB_460K.', '', publication.traits)
publication.traits <- gsub('cov_', '', publication.traits)
publication.traits <- gsub('repro_', '', publication.traits)

rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('UKB_460K.', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('cov_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('repro_', '', rownames(significance.matrix))

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

column.split = c(
    rep('Excitatory', 15),
    rep('Inhibitory', 18),
    rep('Others', 7)
    )

# make color function
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

# select some traits for publication to visualize
publication.traits <- c(
    'Schizophrenia_Pardinas2018',
    'MDD_Howard2019',
    'BIP_Mullins2021',
    'EDU_YEARS',
    'SWB',
    'Alzheimers_Jansen2019'
    )

# create heatmap:
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
    # row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-CpG-cell-type-brain-traits-proportion.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();
