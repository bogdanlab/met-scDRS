### gene-ontology-plots.R #########################################################################
# purpose: plot the gene ontologies for main plots:

### PREAMBLE ######################################################################################
# load in packages:
library(Seurat);
library(sceasy);
library(ggplot2);
library(data.table);
library(ggpubr);
library(clusterProfiler);
library(enrichplot);
library(cowplot);
library(grid);

system.date <- Sys.Date()

# load in the saved gene ontology result:
gene.ontology.result <- readRDS("/u/project/pasaniuc/lixinzhe/R_saves/2024-07-08-75-traits-GO-result.rds")

### VISUALIZATION #################################################################################
# visualize the gene ontology for MDD:
dot.plot <- dotplot(
    gene.ontology.result[['PASS_MDD_Howard2019']],
    showCategory=10
    ) + 
    ggtitle(paste0("MDD gene ontology enrichments")) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(text = element_text(size = 20))


# plot out the umap:
dot.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
plot.name <- paste0(dot.path, system.date, '-MDD-gene-ontology-enrichment-dotplot.png')
png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );
print(dot.plot);
dev.off();

# also make a network plot:
