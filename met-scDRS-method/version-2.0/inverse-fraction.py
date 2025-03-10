### inverse-fraction.py ###########################################################################
# purpose: process the fraction in python, optimized for memory

### PREAMBLE ######################################################################################
# load in libraries:
import argparse
import pandas as pd
import time

# load in user input:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process methylation fraction')
    parser.add_argument('fraction_path', type=str, help='Path to the input CSV file')
    parser.add_argument('output_path', type=str, help='Path to the output processed CSV file')
    args = parser.parse_args()

fraction_path = args.fraction_path.strip()
output_path = args.output_path.strip()

# for script testing:
# fraction_path = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/unique_mch_gene_fraction.csv"
# output_path = "/u/scratch/l/lixinzhe/tmp-file/tmp-output/processed-mch.csv"

# define a custom function to check the data is 1) numeric, and 2) bound between 0 and 1:
def fraction_check(pandas_df):
    # Check if all entries are numeric
    numeric_check = pandas_df.applymap(lambda x: isinstance(x, (int, float))).all().all()
    
    # Check if all entries are between 0 and 1
    bounded_check = pandas_df.applymap(lambda x: 0 <= x <= 1).all().all()
    
    # Output warning if the result is not bounded and not numeric:
    if not numeric_check or not bounded_check:
        raise ValueError("Error: non numeric or not bounded between 0 and 1, check data!")
    return numeric_check & bounded_check

### COMPUTE #######################################################################################
# load in fraction as dask object:
print('loading fraction: ')
fraction = pd.read_csv(fraction_path, index_col = 0)

# perform check for the matrix:
print('performing fraction check:')
print(fraction.head())
fraction_check(fraction)

# compute the variance for all columns
fraction_variance = fraction.var(axis = 0)

# filter based on low variance:
low_variance_quantile = fraction_variance.quantile(0.05)
low_variance_index = fraction_variance < low_variance_quantile

# invert fraction:
fraction *= -1
fraction += 1

# Filter out columns with low variance
processed_fraction = fraction.loc[:, ~low_variance_index]

# check again before output data:
fraction_check(processed_fraction)
print('Finished processing fraction, saving processed fraction')

### OUTPUT ########################################################################################
# put the column name "cell" back:
processed_fraction.insert(0, 'cell', processed_fraction.index)

# output the data:
processed_fraction.to_csv(output_path, index = False)
