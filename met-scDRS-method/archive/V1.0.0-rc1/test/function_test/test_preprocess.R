library(data.table)

python_preprocessed = data.frame(fread('/u/scratch/l/lixinzhe/revision_scratch/v1.0.0-rc1/{today}_inversed_clipped_test.csv'), row.names = 1)

previous_preprocessed = data.frame(fread('/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv'), row.names = 1)

print(dim(previous_preprocessed) == dim(python_preprocessed)) # return true
print(all(colnames(previous_preprocessed) == colnames(python_preprocessed))) # all true
print(all(rownames(previous_preprocessed) == rownames(python_preprocessed))) # true

max_abs_diff <- max(abs(previous_preprocessed - python_preprocessed))
cat('maximum absolute difference between two preprocessing :', max_abs_diff, '\n')
# maximum absolute difference between two preprocessing : 7.327327e-08 

