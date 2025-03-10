### brain-trait-heatmap.R #########################################################################
# purpose: draw a heatmap that looks at the proportion of significant cells per cell type

### PREAMBLE ######################################################################################
# load packages:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# load in data:
significance.matrix <- read.table(
    file = '/u/home/l/lixinzhe/project-geschwind/plot/2024-04-03-mch-MajorType-gse132489-fraction-mch-74traits-cell-class-significance-proportion.csv',
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
cell.type.order <- c(
    all.region,
    ANP.region,
    cortex.region,
    OLF.exc.region,
    HIP.region,
    cortex.hip.olf.region,
    CNU.region,
    OLF.region
    );

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

# perform checks:
stopifnot(length(cell.type.order) == ncol(significance.matrix));
stopifnot(cell.type.order %in% colnames(significance.matrix))

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
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-mouse-cell-type-brain-traits-proportion.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();
