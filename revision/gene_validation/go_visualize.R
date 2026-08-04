# purpose: compute gene ontology and visuazlie the go terms

###########################################################################################
######                                 load in packages                              ######
###########################################################################################
library(clusterProfiler);
library(enrichplot);
library(cowplot);
library(grid);
library(dplyr)
library(ggplot2);


# load in function;
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/'
source(paste0(function.path, 'gene-ontology-caller.R'))
system.date <- Sys.Date();


# define number of top selected genes for computing the go enrichments:
ontology_gene_num <- 100;

# load in the data:
all_cor = read.table("/u/home/l/lixinzhe/project-geschwind/port/scDRS/test_gene_methylation_score_correlations.tsv", sep = '\t', header = TRUE)
gsea.c5 <- read.gmt('/u/project/geschwind/lixinzhe/data/c5.all.v2023.1.Hs.entrez.gmt');

###########################################################################################
######                               Output the Bipolar and MDD Correlation          ######
###########################################################################################
diseases = c('PASS_BIP_Mullins2021', 'PASS_MDD_Howard2019')
for (disease in diseases){
    disease_cor = all_cor[all_cor$trait == disease, ]
    disease_cor = disease_cor[order(disease_cor$corr), ]
    top_genes = tail(disease_cor, ontology_gene_num)
    
    # output:
    file_pth = '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/top_genes_correlation/'
    write.table(top_genes, sep = '\t', row.names = FALSE, file = paste0(file_pth, disease, 'top_correlated_genes.tsv'))
}

###########################################################################################
######                                    Make Plots                                 ######
###########################################################################################

# for the highest correlation genes, get the gene ontology:
diseases = unique(all_cor$trait)
for (disease in diseases){
    # get the disease correlation:
    disease_cor = all_cor[all_cor$trait == disease, ]
    disease_cor = disease_cor[order(disease_cor$corr), ]
    
    # get the disease correlation:
    top_genes = tail(disease_cor, ontology_gene_num)$gene
    background_genes = unique(disease_cor$gene)
    
    # call the ontology:
    gene.ontology.result <- gene.ontology.caller(
        x = top_genes,
        background = background_genes,
        terms = gsea.c5,
        visualize = TRUE
        );
    readable.result <- setReadable(gene.ontology.result, 'org.Hs.eg.db', 'ENTREZID')
    
    # visualize the gene ontology for MDD:
    dot.plot <- dotplot(
        gene.ontology.result,
        showCategory=10
        ) + 
        ggtitle(paste0("gene ontology for ", disease)) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(text = element_text(size = 20))
    
    # plot out the umap:
    dot.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/dot/'
    plot.name <- paste0(dot.path, system.date, '-', disease, '-gene-ontology-enrichment-dotplot.png')
    png(
        filename = plot.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    print(dot.plot);
    dev.off();
    
    # deal with the network:
    cor_value = tail(disease_cor, ontology_gene_num)$corr
    names(cor_value) = tail(disease_cor, ontology_gene_num)$gene
    network.plot <- cnetplot(
        readable.result,
        categorySize = "pvalue",
        color.params = list(foldChange = cor_value, edge = TRUE),
        circular = TRUE) + 
        ggtitle(disease) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(legend.text = element_text(size = 12))
        # theme(legend.position = "none")
    
    network.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/network/'
    plot.name <- paste0(network.path, system.date, '-', disease, '-gene-ontology-enrichment-network-legend.png')
    png(
        filename = plot.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    legend <- cowplot::get_legend(network.plot);
    grid.newpage()
    grid.draw(legend)
    dev.off()
    
    # also get the network itself:
    plot.name <- paste0(network.path, system.date, '-', disease, '-gene-ontology-enrichment-network.png')
        network.plot <- network.plot + theme(legend.position = "none")
        png(
            filename = plot.name,
            width = 10,
            height = 10,
            units = 'in',
            res = 400
            );
        print(network.plot);
        dev.off();
}
    
    
    