### visualization.R ###############################################################################
# purpose: visualize the results on scDRS run using AIBS

### PREAMBLE ######################################################################################
# load in libraries:
library(data.table)
require(ggplot2);
require(ComplexHeatmap);
require(circlize);
system.date <- Sys.Date()

# load in paths:
scDRS_directory <- '/u/project/geschwind/lixinzhe/scDRS-output/scDRS-output/AIBS-psych-trait-scDRS/with_cov/'

# load in scDRS for MDD:
mdd <- fread(paste0(scDRS_directory, 'PASS_MDD_Howard2019.score.gz'), sep = '\t', data.table = FALSE)
mdd <- data.frame(mdd, row.names = 1)

# load in the meta file:
meta_file <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/metadata.csv'
meta <- data.frame(fread(meta_file, sep = ',', data.table = FALSE), row.names = 1)
meta <- meta[rownames(mdd), ]

### PROCESS #######################################################################################
# crate a class label:
meta$class <- 'Non-neuronal'
meta[meta$class_label == 'Glutamatergic', 'class'] = 'Excitatory'
meta[meta$class_label == 'GABAergic', 'class'] = 'Inhibitory'

# find the set of cells that are marginally significant:
marginal_significant <- rownames(mdd)[mdd$pval < 0.05]

# find the identify of the set of cells that are marginally significant
table(meta[marginal_significant, 'subclass_label'], meta[marginal_significant, 'cortical_layer_label'])

# find the identify of the set of cells that are marginally significant
table(meta[marginal_significant, 'class_label'], meta[marginal_significant, 'cortical_layer_label'])

# make heatmap that shows the marginal significant cell count:
# create a color function, note 100 is set as incremental color
significance.matrix <- table(meta[marginal_significant, 'class'], meta[marginal_significant, 'cortical_layer_label'])
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
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-marginal-sig-count.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 30, 110), "mm"));
dev.off();

# also make a bar chart on the significance cell count in different class of neurons:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-AIBS-marginal-sig-count-class.png')
png(
    filename = output.path,
    width = 4,
    height = 4,
    units = 'in',
    res = 400
    );

# make a bar plot:
freq_df <- as.data.frame(table(meta[marginal_significant, 'class']))
ggplot(freq_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "AIBS atlas dataset", x = "cell class", y = "count") +
  theme_classic()

dev.off()