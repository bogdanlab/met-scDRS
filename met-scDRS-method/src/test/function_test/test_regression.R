library(data.table)

python_preprocessed = data.frame(fread('/u/scratch/l/lixinzhe/revision_scratch/batch_regress/scdrs_regress_fx.csv'), row.names = 1)

previous_preprocessed = data.frame(fread('/u/scratch/l/lixinzhe/revision_scratch/batch_regress/scdrs_regress_inplace_fx.csv'), row.names = 1)

print(dim(previous_preprocessed) == dim(python_preprocessed)) # return TRUE TRUE
print(all(colnames(previous_preprocessed) == colnames(python_preprocessed))) # return TRUE
print(all(rownames(previous_preprocessed) == rownames(python_preprocessed))) # return TRUE

max_abs_diff <- max(abs(previous_preprocessed - python_preprocessed))
cat('maximum absolute difference between two preprocessing :', max_abs_diff, '\n')
# maximum absolute difference between two preprocessing : 0.00013858


