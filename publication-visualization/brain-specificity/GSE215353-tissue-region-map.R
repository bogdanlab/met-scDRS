### GSE215353-tissue-region-map.R #################################################################
# purpose: look at the contingency table between region (verbal) and AIBS code

### PREAMBLE ######################################################################################
library(ComplexHeatmap)
system.date <- Sys.Date()

# load in data:
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
output.path <- "/u/home/l/lixinzhe/project-geschwind/plot/"
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

### ANALYSIS ######################################################################################
# make contingency table:
wide.plot.df <- table(meta$tissue, meta$X_Region)

plot <- Heatmap(
    wide.plot.df,
    name = 'cell count',
    rect_gp = gpar(col = "black", lwd = 2),
    width = unit(10 * ncol(wide.plot.df),"mm"),
    height = unit(10 * nrow(wide.plot.df),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    cluster_rows = FALSE,
    cluster_columns = FALSE
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
plot.path <- paste0(output.path, system.date, '-GSE215353-tissue-region-count-heatmap.png')
png(
    filename = plot.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();
