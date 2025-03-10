### score properties analysis.R ################################################################
# purpose: computes the average of the met-scDRS for regions within L2/3 MDD significant cells
# computes the average number of significant cells in brain vs non brain traits

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)

# specify paths:
scDRS.directory <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/';
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv';
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date()

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# read trait info:
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# load in data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# create progress bar:
cat('loading scores \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

# read into the empty list:
for (result in score.files) {
    risk.score[[result]] <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
        );
    risk.score[[result]] <- data.frame(risk.score[[result]], row.names = 1)
    pb$tick()
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

### SUMMARIZATION MDD #############################################################################
# for MDD, within L2/3 neurons, summarize the average score by region
mdd.score <- risk.score[["PASS_MDD_Howard2019"]];
mdd.score$region <- meta$tissue[match(rownames(mdd.score), rownames(meta))]
mdd.score$cell_type <- meta$X_MajorType[match(rownames(mdd.score), rownames(meta))]
l23.mdd.score <- mdd.score[mdd.score$cell_type == "L2/3-IT", ]

# summarize by average and sd of each region
summarize.table <- l23.mdd.score %>% group_by(region) %>% summarize(average = mean(norm_score), sd = sd(norm_score))
summarize.table <- data.frame(summarize.table)

# output this summarize table as a supplementary table:
write.table(
    x = summarize.table,
    file = paste0(
        '/u/home/l/lixinzhe/project-pasaniuc/plot/',
        system.date,
        '-l23-MDD-region-summary.csv'
        ),
    sep = ',',
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
    )

### SUMMARIZE TRAITS ##############################################################################
# Find the average number of significant cells that are brain vs non brain traits:
significant.number <- data.frame(
    trait = names(risk.score),
    category = trait.info$Category[match(names(risk.score), trait.info$Trait_Identifier)],
    number = NA
    )
rownames(significant.number) <- significant.number$trait
for (disease in names(risk.score)) {
    fdr.p <- p.adjust(risk.score[[disease]]$pval, method = 'fdr')
    significant.number[disease, 'number'] <- sum(fdr.p < 0.1)
    }

# output the significant number table:
write.table(
    x = significant.number,
    file = paste0(
        '/u/home/l/lixinzhe/project-pasaniuc/plot/',
        system.date,
        '-number-of-significant-cells-by-trait-summary.csv'
        ),
    sep = ',',
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
    )

# summarize by trait categories:
summarize.table <- data.frame(significant.number %>% group_by(category) %>% summarize(average = mean(number), sd = sd(number)))

write.table(
    x = summarize.table,
    file = paste0(
        '/u/home/l/lixinzhe/project-pasaniuc/plot/',
        system.date,
        '-summarize-significant-cells-by-trait-categories.csv'
        ),
    sep = ',',
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
    )
