"""
### DEVELOPING ###
Updated on July 19th, 23

Description:
This file is for wrangling the sourmash output to create the count table
"""

# Import packages:
import glob
import os
import argparse
import pandas as pd

#############
# ARGUMENTS #
#############

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str)

args = parser.parse_args()
PATH = args.input

##############
# PARAMETERS #
##############

SCALE = 1000  # scale number of the sourmash sketch

##################
# DATA WRANGLING #
##################

# Concatenate all lineages.csv files:
joined_files = os.path.join(PATH, "*.csv")
joined_list = glob.glob(joined_files)
li = []


for filename in joined_list:
    alt_df = pd.read_csv(filename,
                         usecols=['query_name',
                                  'name',
                                  'intersect_bp',
                                  'unique_intersect_bp',
                                  'average_abund',
                                  'match_containment_ani'])
    # suggested filter by Zhenfeng
    fil_df = alt_df[(alt_df['match_containment_ani'] >= 0.935) |
                    (alt_df['intersect_bp'] >= SCALE*1000)]

    minor_df = fil_df.drop(columns=['match_containment_ani', 'intersect_bp'])
    li.append(minor_df)

df = pd.concat(li, axis=0, ignore_index=True)

# Compute the Ident counts in each sample:
df['n_unique_kmer'] = (df['unique_intersect_bp']/SCALE)*df['average_abund']
df = df[['query_name', 'name', 'n_unique_kmer']]
df[['Ident', 'Ident_name']] = df["name"].str.split(" ", n=1, expand=True)

########################################################################
# FINAL STEP: Extract the 'Ident' counts, and pivot into the ASV table:#
########################################################################

# Create the Count table:
df2 = pd.pivot_table(df,
                     index='query_name',
                     columns='Ident',
                     values='n_unique_kmer',
                     fill_value=0)

df2 = df2.astype('int')

print('#' * len(f'# The count table is created by {len(df2)} samples #'))
print(f'# The count table is created by {len(df2)} samples #')
print('#' * len(f'# The count table is created by {len(df2)} samples #'))
print(df2.head())

df2.to_csv('/data/namlhs/gmrepo/data/alternative/absolute_count_raw.csv',
           sep='\t')
