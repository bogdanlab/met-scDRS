### cell-type-tissue-mc-model-distribution.R ######################################################
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

# Extract the arguments
args <- commandArgs(trailingOnly = TRUE)
model.result.dir <- args[1]

# load in the data:
model.result.dir <- "/u/scratch/l/lixinzhe/tmp-file/region-cell-type-analysis-sig-only/PASS_MDD_Howard2019/"
mc.results <- list.files(model.result.dir, pattern = '*-model-summary.csv');

# create a progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(mc.results),
    clear = FALSE
    );

cat('reading in result files \n')
result.collection <- NULL;
for (result.file in mc.results) {
    pair.index <- gsub('-MC-linear-model-summary.csv', '', result.file)
    result.collection[[pair.index]] <- data.frame(
            fread(
            file = paste0(model.result.dir, result.file),
            sep = ',',
            data.table = FALSE
            ),
        row.names = 1
        )
    pb$tick()
    }

### COMPUTE PERMUTED VALUE ########################################################################
# compute the permuted value:
cat('compute permuted value: \n')

# create a placeholder:
permuted.value <- rep(NA, length = length(result.collection))
names(permuted.value) <- names(result.collection)
beta <- zscore <- permuted.value

# compute the p value:
for (result in names(result.collection)) {
    foreground.effect <- result.collection[[result]]['foreground', 'beta']
    if (is.na(foreground.effect)) {
        permuted.value[result] <- NA;
        beta[result] <- NA;
        } else {
            numerator <- 1 + sum(abs(result.collection[[result]][, 'beta']) > abs(foreground.effect))
            denominator <- nrow(result.collection[[result]])
            pval <- numerator / denominator
            permuted.value[result] <- pval
            beta[result] <- result.collection[[result]]["foreground", 'beta']

            # also compute the z score using (foreground - background mean) / background sd
            # grep out the control mc draws:
            control.index <- grep('ctrl_norm_score_*', rownames(result.collection[[result]]))
            z.mean <- mean(result.collection[[result]][control.index, 'beta'])
            z.sd <- sd(result.collection[[result]][control.index, 'beta'])
            zscore[result] <- (foreground.effect - z.mean) / z.sd
        }
    }

### VISUALIZATOIN #################################################################################
# visualize the data:
permuted.value.complete <- permuted.value[!is.na(permuted.value)]
beta.complete <- beta[!is.na(beta)]
neg.logp <- -log10(permuted.value.complete)
z.complete <- zscore[!is.na(zscore)]

# construct plotting data frame:
plot.df.complete <- data.frame(
    neg_log_p = neg.logp,
    effect = beta.complete,
    zscore = z.complete,
    cell_type_region_pairs = names(neg.logp),
    cell_type = gsub(':.*', '', names(neg.logp)),
    region = gsub('.*:', '', names(neg.logp))
    );

plot.df <- data.frame(
    neg_log_p = -log10(permuted.value),
    effect = beta,
    zscore = zscore,
    cell_type_region_pairs = names(beta),
    cell_type = gsub(':.*', '', names(beta)),
    region = gsub('.*:', '', names(beta))
    )

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
color.dictionary <- data.frame(cell.type, cell.class)
plot.df$cell_class <- color.dictionary$cell.class[match(plot.df$cell_type, color.dictionary$cell.type)]

### VISUALIZE #####################################################################################
# convert data format from long to wide:
wide.plot.df <- reshape(
    plot.df[, c('region', 'cell_type', 'zscore')],
    timevar = 'cell_type',
    idvar = 'region',
    direction = 'wide'
    )
row.names(wide.plot.df) <- wide.plot.df$region
wide.plot.df <- wide.plot.df[, -1]
colnames(wide.plot.df) <- gsub('zscore.', '', colnames(wide.plot.df))

# make a heatmap:
col.fun <- colorRamp2(
    c(
        min(plot.df$zscore, na.rm = TRUE),
        0,
        max(plot.df$zscore, na.rm = TRUE)
        ),
    c('#0571b0', 'white', '#ca0020')
    );
heatmap.legend.param <- list(
    at = c(
        round(
            min(plot.df$zscore, na.rm = TRUE),
            digit = 2 # round to 2 digit
            ),
        0,
        round(
            max(plot.df$zscore, na.rm = TRUE),
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
wide.plot.df <- wide.plot.df[, cell.type]
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
    "L6-IT"
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

# Now create the heatmap:
col.fun <- colorRamp2(
    c(
        min(wide.plot.df, na.rm = TRUE),
        0,
        max(wide.plot.df, na.rm = TRUE)
        ),
    c('#0571b0', 'white', '#ca0020')
    );
selected.plot.df <- wide.plot.df[selected.region, selected.cell.type]
plot <- Heatmap(
    as.matrix(selected.plot.df),
    name = 'zscore',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = rownames(selected.plot.df),
    column_order = colnames(selected.plot.df),
    width = unit(10 * ncol(selected.plot.df),"mm"),
    height = unit(10 * nrow(selected.plot.df),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    cluster_rows = FALSE,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

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
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
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
