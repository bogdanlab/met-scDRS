### common-axis-cell-type-region-predictiveness-heatmap.R #########################################
# purpose: calculate the permuted p value for the cell type tissue mc model distribution
# this script will visualize the cell type - tissue pairs for one disease

### PREAMBLE ######################################################################################
# load libraries:
require(docopt)
require(circlize);
require(ComplexHeatmap);
require(ggplot2);
require(tidyverse);
require(ggrepel);
require(data.table);

# define parameters:
output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date();

# load in data:
# load in meta data:
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# load in the scDRS score:
scDRS.file <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.score.gz"
# load in the data
full.score <- fread(
    file = scDRS.file,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE
    );
rownames(full.score) <- full.score$cell

### ANALYSIS ######################################################################################
# get the set of cells that we wish to plot:
all(rownames(meta) == rownames(full.score))
meta$score = full.score$zscore
cell_type_tissue_average = meta %>% 
    group_by(X_MajorType, tissue) %>% 
    summarize(average = mean(score), sd = sd(score), cell_num = length(score)) %>% data.frame()
filtered = cell_type_tissue_average[cell_type_tissue_average$cell_num > 10, ]

# pivot to wider format
heatmap_df = filtered[, c('X_MajorType', 'tissue', 'average')] %>% pivot_wider(names_from = X_MajorType, values_from = average) %>% data.frame(row.names = 1, check.names = FALSE)

# create a dictionary for plotting the color on excitatory, inhibitory and non neuronal:
cell.type <- c(
    "L2-3-IT",
    "L4-IT",
    "L5-ET",
    "L5-IT",
    "L5-6-NP",
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
cell.class <- c(
    rep('Excitatory', 15),
    rep('Inhibitory', 18),
    rep('Others', 7)
    )
colnames(heatmap_df) = gsub('/', '-', colnames(heatmap_df))
color.dictionary <- data.frame(cell.type, cell.class)
cell_class <- color.dictionary$cell.class[match(colnames(heatmap_df), color.dictionary$cell.type)]

### VISUALIZE #####################################################################################
# make a heatmap:
col.fun <- colorRamp2(
    c(
        min(filtered$average, na.rm = TRUE),
        0,
        max(filtered$average, na.rm = TRUE)
        ),
    c('#0571b0', 'white', '#ca0020')
    );
heatmap.legend.param <- list(
    at = c(
        round(
            min(filtered$average, na.rm = TRUE),
            digit = 2 # round to 2 digit
            ),
        0,
        round(
            max(filtered$average, na.rm = TRUE),
            digit = 2 # round to 2 digit
            )
        )
    );
column.split = c(
    rep('Excitatory', 15),
    rep('Inhibitory', 18),
    rep('Others', 7)
    )
    
# Now create the heatmap:
wide.plot.df <- heatmap_df[, cell.type]
plot <- Heatmap(
    as.matrix(wide.plot.df),
    name = 'zscore',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    # row_order = publication.traits,
    column_order = colnames(wide.plot.df),
    width = unit(10 * ncol(wide.plot.df),"mm"),
    height = unit(10 * nrow(wide.plot.df),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    cluster_rows = FALSE,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
plot.path <- paste0(output.path, system.date, '-region-cell-type-zscore-heatmap.png')
png(
    filename = plot.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

write.table(
    wide.plot.df,
    file = gsub('png', 'csv', plot.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

### VISUALIZE SELECTED CELL TYPE AND REGION #######################################################
selected.cell.type <- c(
    "L2-3-IT",
    "L4-IT",
    "L5-IT",
    "L6-IT",
    "MGC",
    "OPC"
    )
selected.region <- c(
    "primary-visual-cortex",
    "secondary-visual-cortex",
    "primary-somatosensory-cortex",
    "cuneus-cortex",
    "primary-motor-cortex",
    "parietal-cortex",
    "primary-auditory-cortex",
    "Brodmann-(1909)-area-19",
    "prefrontal-cortex",
    "Brodmann-(1909)-area-46",
    "granular-insular-cortex",
    "inferior-temporal-gyrus",
    "middle-temporal-gyrus",
    "posterior-parahippocampal-gyrus",
    "Brodmann-(1909)-area-25",
    "insula",
    "anterior-cingulate-cortex",
    "Brodmann-(1909)-area-38",
    "agranular-insular-cortex",
    "piriform-cortex",
    "medial-entorhinal-cortex",
    "lateral-entorhinal-cortex"
    )
selected.region = gsub('-', ' ', selected.region)
selected.plot.df <- wide.plot.df[selected.region, selected.cell.type]

# Now create the heatmap:
col.fun <- colorRamp2(
    c(
        min(filtered$average, na.rm = TRUE),
        0,
        max(filtered$average, na.rm = TRUE)
        ),
    c('#0571b0', 'white', '#ca0020')
    );
heatmap.legend.param <- list(
    at = c(
        round(
            min(filtered$average, na.rm = TRUE),
            digit = 2 # round to 2 digit
            ),
        0,
        round(
            max(filtered$average, na.rm = TRUE),
            digit = 2 # round to 2 digit
            )
        )
    );

plot <- Heatmap(
    as.matrix(selected.plot.df),
    name = 'average\nmet-scdrs',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = rownames(selected.plot.df),
    column_order = colnames(selected.plot.df),
    width = unit(10 * ncol(selected.plot.df),"mm"),
    height = unit(10 * nrow(selected.plot.df),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    cluster_rows = FALSE,
    heatmap_legend_param = heatmap.legend.param,
    show_heatmap_legend = FALSE
    );
plot.size <- draw(plot, padding = unit(c(10, 10, 10, 30), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
plot.path <- paste0(output.path, system.date, '-region-cell-type-zscore-selected-heatmap.png')
png(
    filename = plot.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, padding = unit(c(10, 10, 10, 30), "mm"));
dev.off();

# next draw the legend:
lgd = Legend(
    col_fun = col.fun, 
    at = c(
        round(
            min(selected.plot.df, na.rm = TRUE),
            digit = 2 # round to 2 digit
            ),
        0,
        round(
            max(selected.plot.df, na.rm = TRUE),
            digit = 2 # round to 2 digit
            )
        ),
    title = 'average\nmet-scDRS',
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
    )
plot.path <- paste0(output.path, system.date, '-region-cell-type-zscore-selected-heatmap-legend.png')
png(
    filename = plot.path,
    width = 3,
    height = 3,
    units = 'in',
    res = 400
    );
draw(lgd);
dev.off();

# write out the selected matrix:
write.table(
    selected.plot.df,
    file = gsub('png', 'csv', plot.path),
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

# load in the file and calculate the average of average:
selected.plot.df = read.table(file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-01-09-region-cell-type-zscore-selected-heatmap.csv', sep =',')
mean(as.matrix(selected.plot.df)[c('prefrontal-cortex','Brodmann-(1909)-area-46', 'inferior-temporal-gyrus','middle-temporal-gyrus'),], na.rm = T) # 1.072
sd(as.matrix(selected.plot.df)[c('prefrontal-cortex','Brodmann-(1909)-area-46', 'inferior-temporal-gyrus','middle-temporal-gyrus'),], na.rm = T) # 0.673

full.plot.df = read.table(file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-01-09-region-cell-type-zscore-heatmap.csv', sep = ',')
mean(as.matrix(full.plot.df)[c('prefrontal-cortex','Brodmann-(1909)-area-46', 'inferior-temporal-gyrus','middle-temporal-gyrus'), grep('IT', colnames(full.plot.df))], na.rm = T) # 1.072
