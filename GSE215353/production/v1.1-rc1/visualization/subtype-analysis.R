### subtype-analysis.R ############################################################################
# purpose: investigate how our data varies in subtype categories:

### PREAMBLE ######################################################################################
# load in libraries:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);
library(Seurat);

# load in data:
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

subtype_meta <- data.table::fread(
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/table-S5-cell-meta.csv',
    sep = ',',
    )
subtype_meta <- data.frame(subtype_meta, row.names = 1)

# read in the results:
trait.score <- read.table(
    file = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.score.gz',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    )
trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr')
system.date <- Sys.Date()
p.cutoff <- 0.1

# load in all the traits z scores:
# read data:
score.files <- list.files('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/', pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

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
list.names <- gsub('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/', '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

### PROCESS #######################################################################################
# merge the subtype information onto the actual meta:
meta$subtype <- subtype_meta[rownames(meta), 'SubType']

# make a umap on the L2/3 neurons colored by subtype:
# plot out the spatial region for the L2/3 neurons:
for (cell_type in unique(meta$X_MajorType)){
    cat('processing', cell_type, '\n')
    cell_type_cells <- rownames(meta)[meta$X_MajorType == cell_type];
    cell_type_meta <- meta[cell_type_cells, ]
    cell_type_meta$subtype <- paste0(cell_type, '_', as.numeric(as.factor(cell_type_meta$subtype)))
    
    if (length(unique(cell_type_meta$subtype)) > 1){  
        plot.df <- data.frame(
            UMAP_1 = cell_type_meta$UMAP_1,
            UMAP_2 = cell_type_meta$UMAP_2,
            tissue = as.factor(cell_type_meta$tissue),
            cell_type = as.factor(cell_type_meta$X_MajorType),
            met_scdrs = trait.score[cell_type_cells, 'zscore'],
            fdr = trait.score[cell_type_cells, 'fdr'],
            sub_type = as.factor(cell_type_meta$subtype)
            );
        rownames(plot.df) = rownames(cell_type_meta);

        # Create plot:
        gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = sub_type)) +
            geom_point() +
            theme_classic() +
            ggtitle(paste0('GSE215353 ', cell_type, ' subtypes UMAP')) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab('UMAP1') +
            ylab('UMAP2') +
            theme(legend.position="none") +
            theme(text = element_text(size = 20))
        # Create plot:
        gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = sub_type)) +
            geom_point() +
            theme_classic() +
            ggtitle(paste0('GSE215353 ', cell_type, ' subtypes UMAP')) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab('UMAP1') +
            ylab('UMAP2') +
            theme(legend.position="none") +
            theme(text = element_text(size = 20))

        gplot.label <- LabelClusters(plot = gplot, id = 'sub_type', col = 'black', size = 5)
        plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/subtype-analysis/", system.date, '-GSE215353-', gsub('/','-', cell_type), '-subtype-with-legend-umap.png')
        png(
            filename = plot.path,
            width = 14,
            height = 14,
            units = 'in',
            res = 400
            );
        print(gplot.label)
        dev.off();

        # also draw out the met-scDRS:
        significant.cell <- rownames(plot.df)[plot.df$fdr < p.cutoff]
        insignificant.cell <- setdiff(rownames(plot.df), significant.cell)

        gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2)) +
            geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
            geom_point(data = plot.df[significant.cell, ], aes(colour = met_scdrs)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(paste0(cell_type, ' MDD Howard 2019')) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab('UMAP 1') +
            ylab('UMAP 2') +
            labs(color = "met-scDRS") +
            theme(text = element_text(size = 20))

        # draw out the plot:
        output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/subtype-analysis/'
        output.path <- paste0(output.dir, system.date, '-MDD-scDRS-score-umap-', gsub('/','-', cell_type), '-only.png')
        png(
            filename = output.path,
            width = 14,
            height = 14,
            units = 'in',
            res = 400
            );
        print(gplot)
        dev.off();

        ### VISUALIZATION - proportion heatmap ############################################################
        # for each of the result, find out the number of significant cells in each cell type:
        group.index <- 'subtype'
        cell.type <- unique(cell_type_meta[, group.index]);
        significance.matrix <- matrix(NA, nrow = length(risk.score), ncol = length(cell.type));
        rownames(significance.matrix) <- names(risk.score);
        colnames(significance.matrix) <- cell.type;
        trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

        # create the traits by cell type matrix:
        plot.type <- 'proportion'
        for (result in names(risk.score)) {
            # find the index for significant cells:
            significant.cell <- rownames(risk.score[[result]])[
                p.adjust(risk.score[[result]]$pval, method = 'fdr') < p.cutoff
                ];

            # for each cell type, check number of significant cells are part of that cell type:
            for (type in cell.type) {
                # locate cell id that belong in the cell type:
                cell.type.cell <- rownames(cell_type_meta)[cell_type_meta[, group.index] %in% type];

                # find the number of cells that are in each of the cell type category:
                if (plot.type == 'count') {
                    significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell);
                    } else {
                        significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell) /
                            length(cell.type.cell);
                    }
                }
            }

        # draw out the proportion heatmap:
        # select traits to plot:
        publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

        # reformat traits:
        publication.traits <- gsub('PASS_', '', publication.traits)
        rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
        cell.type.order <- paste0(cell_type, '_', seq(1, ncol(significance.matrix)))

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
            heatmap_legend_param = heatmap.legend.param
            );
        plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

        # measure the size of the heatmap:
        heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
        heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

        # use the measured width and height for drawing:
        output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/subtype-analysis/', system.date, '-', gsub('/','-', cell_type), '-subtype-brain-traits-proportion.png')
        png(
            filename = output.path,
            width = heatmap.width,
            height = heatmap.height,
            units = 'in',
            res = 400
            );
        draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
        dev.off();
        
        # also generate the supplementary tables:
        table.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/subtype-analysis/', system.date, '-', gsub('/','-', cell_type), '-subtype-brain-traits-proportion.tsv')
        write.table(significance.matrix[publication.traits, cell.type.order], file=table.path, sep = '\t')
        }
    }
   

### chisquare test between subtype and region #####################################################
# for all subtype:
contingency.table <- table(cell_type_meta$subtype, cell_type_meta$tissue)
print(chisq.test(contingency.table, simulate.p.value = TRUE, B = 200000))
