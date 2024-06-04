"""
### DEVELOPING ###
Updated on June 19, 23

Description:
This file creates the Tax Table (ASV dictionary in specific taxonomic levels)
"""

# Import packages
import pandas as pd
from IPython.display import display

# Import data from SOURMASH references
tax_data = pd.read_csv("/data/namlhs/gmrepo"
                       "/gtdb-rs214.lineages.csv",
                       low_memory=False)

# tax_table = tax_data.drop(columns=['taxid', 'strain'])
tax_data.columns = tax_data.columns.str.capitalize()
tax_data.rename(columns={'Superkingdom': 'Domain'}, inplace=True)

# Assembly the tax table into the asv
tax_data['ASV_cluster'] = tax_data[['Domain',
                                    'Class',
                                    'Order',
                                    'Family',
                                    'Genus',
                                    'Species']].apply(';'.join,
                                                      axis=1)

asv_counts = tax_data.groupby('ASV_cluster').size().reset_index(name='Count')
asv_counts['ASV'] = ['ASV'+str(i) for i in range(asv_counts.shape[0])]

asv_dict = dict(zip(asv_counts.ASV_cluster, asv_counts.ASV))

tax_data['ASV'] = tax_data['ASV_cluster'].map(asv_dict)
tax_data = tax_data.drop(columns=['ASV_cluster'])

display(tax_data)

tax_data.to_csv('/data/namlhs/nextflow-sourmash/tax_table.csv',
                sep='\t')
