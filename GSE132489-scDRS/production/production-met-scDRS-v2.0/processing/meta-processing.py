### meta-processing.py ############################################################################
# purpose: read in the excel file and output the processed excel:

# load in libraries:
import pandas

# set global options:
pandas.set_option('display.max_columns', None)

# load in the meta data:
meta_data_path = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
meta = pandas.read_excel(meta_data_path, skiprows = 15, index_col = 0)
# meta.shape returns (110294, 27) in shape

# check how many non QC passing cell is in the meta data:
# meta['Pass QC'].value_counts() 
# return True     103982
# False      6312

### OUTPUT ########################################################################################
# write the data:
output_dir = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/"
output_path = output_dir + "all_meta_data.csv"
meta.insert(0, 'cell', meta.index)
meta.to_csv(output_path, index=False)
