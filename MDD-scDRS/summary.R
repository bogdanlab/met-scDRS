### summary.R #####################################################################################
# get summary from the scDRS result to write paper

### PREAMBLE ######################################################################################
# load in libraries:
library(data.table)
require(ggplot2);
require(ComplexHeatmap);
require(circlize);
system.date <- Sys.Date()

# load in paths:
scDRS_directory <- '/u/project/geschwind/lixinzhe/scDRS-output/MDD/'

# load in scDRS for MDD:
mdd <- fread(paste0(scDRS_directory, 'PASS_MDD_Howard2019.score.gz'), sep = '\t', data.table = FALSE)
mdd <- data.frame(mdd, row.names = 1)

# load in the meta file:
meta <- read.table(
    sep = '\t',
    header = TRUE,
    row.names = 1,
    file = paste0('/u/project/geschwind/lixinzhe/data/2023-05-07-MDD-GSE144136-metadata.txt')
    );

### PROCESS #######################################################################################
# put neuronal class into the dataframe:
meta$class <- 'Non-neuronal'
meta[meta$orig.ident == 'Ex', 'class'] = 'Excitatory'
meta[meta$orig.ident == 'Inhib', 'class'] = 'Inhibitory'

# marginally significant p value:
marginal_significant <- rownames(mdd)[mdd$pval < 0.05]
significance.matrix <- table(meta[marginal_significant,'class'], meta[marginal_significant,'cell.type'])
col.fun <- colorRamp2(
    c(
        0,
        max(significance.matrix)
        ),
    c('white','#de2d26')
    );

# create heatmap legend:
heatmap.legend.param <- list(
    at = c(
        0,
        round(
            max(significance.matrix),
            digit = -2 # round to nearest hundreds
            )
        )
    );

# make the plot:
plot <- Heatmap(
    as.matrix(significance.matrix),
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    width = unit(10 * ncol(significance.matrix),"mm"),
    height = unit(10 * nrow(significance.matrix),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 30, 110), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-marginal-sig-count-cell-class-nagy-etal.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 30, 110), "mm"));
dev.off();


freq_df <- as.data.frame(table(meta[marginal_significant, 'class']))

# make a plot on number of significant cells in different classes:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE144136-marginal-sig-count-class.png')
png(
    filename = output.path,
    width = 4,
    height = 4,
    units = 'in',
    res = 400
    );

ggplot(freq_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "case/control dataset", x = "cell class", y = "count") +
  theme_classic()
dev.off()

# make a heatmap:

