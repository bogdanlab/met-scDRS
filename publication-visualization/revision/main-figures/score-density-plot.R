library(ggplot2)
library(ggplot2)
library(dplyr)
library(data.table)

# define parameters:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date();

# load in the MDD risk score
# load in the scDRS score:
scDRS.file <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.full_score.gz"
# load in the data
full.score <- fread(
    file = scDRS.file,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE
    );
rownames(full.score) <- full.score$cell

# faltten, select all normalized score
score_cols <- paste0("ctrl_norm_score_", 0:999)

# Flatten the selected columns into one numeric vector
scores <- unlist(full.score[nrow(full.score), score_cols, drop = FALSE], use.names = FALSE)
scores <- as.numeric(scores)
scores <- scores[is.finite(scores)]

gplot = ggplot(data.frame(score = scores), aes(x = score)) +
  geom_density(color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "control score", y = "Density")

# add a pdf output:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-met-scDRS-MDD-density.pdf')
pdf(
    file = output.path,
    width = 5,
    height = 5
    );
print(gplot)
dev.off();
