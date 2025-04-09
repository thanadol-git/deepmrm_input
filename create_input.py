import pandas as pd
from pyteomics import parser
import re

# Read input files

# Sequence file
qReps = pd.read_csv("data/ComplemEdge_seq_HPRR.csv")
# Method file 
vebios_method = pd.read_csv("data/Vebios_ComplemEdge_method_20241108.csv")

# Extract unique sequences from qReps
qReps['aa_seq'] = qReps['aa_seq'].str.replace('[^A-Z]', '', regex=True)
all_seq = qReps['aa_seq'].unique()

# Function to perform in silico digestion
def digest_df(seq):
    dig_seq = parser.cleave(seq, parser.expasy_rules['trypsin'], missed_cleavages=2)
    dig_seq_df = pd.DataFrame(dig_seq, columns=['peptide'])
    dig_seq_df['aa_seq'] = seq
    Prot_df = qReps[qReps['aa_seq'] == seq].merge(dig_seq_df, on='aa_seq', how='left')
    return Prot_df

# Run digestion on all sequences and combine results
all_dig = pd.concat([digest_df(seq) for seq in all_seq]).drop_duplicates(subset=['peptide'])

# Extract unique peptides from vebios_method
vebios_method['Unique_Seq'] = vebios_method['Compound'].str.replace('[^A-Z]', '', regex=True)
method_seq = vebios_method['Unique_Seq'].unique()

print(vebios_method.head())
# Filter peptides not in method_seq
# # Filter peptides not in method_seq
# filtered_peptides = all_dig[~all_dig['peptide'].isin(method_seq)]['peptide'].unique()

# # # Process vebios_method for deepmrm
# vebios_deepmrm = vebios_method[['Compound', 'Precursor (m/z)', 'Product (m/z)']].copy()
# vebios_deepmrm.columns = ['peptide_id', 'precursor_mz', 'product_mz']
# vebios_deepmrm['is_heavy'] = vebios_deepmrm['peptide_id'].str.contains('heavy')

# # Reformat peptide sequences
# vebios_deepmrm['peptide_id'] = vebios_deepmrm['peptide_id'].str.replace(r'\[.*?\]', '', regex=True)
# vebios_deepmrm['peptide_id'] = vebios_deepmrm['peptide_id'].str.replace(r'\(.*?\)', '', regex=True)
# vebios_deepmrm['peptide_id'] = vebios_deepmrm['peptide_id'].str.replace('C', 'C[+57]')
# vebios_deepmrm['peptide_id'] = vebios_deepmrm['peptide_id'].str.replace(r'\s+', '', regex=True)


# # Find odd numbers of peptides
# count rows by peptide_id and is_heavy in vebios_deepmrm

# peptide_counts = vebios_deepmrm.groupby(['peptide_id', 'is_heavy']).size()
# Extract odd counts from peptie_counts
# odd_counts = peptide_counts[(peptide_counts['False'] %2 == 1)].sum(axis=1)
# odd_pept = peptide_counts[(peptide_counts[False] != peptide_counts[True])].index

# # Extract odd peptides
# odd_peptides = vebios_deepmrm[vebios_deepmrm['peptide_id'].isin(odd_pept)].sort_values(by=['peptide_id', 'is_heavy'])

# # Export results
# vebios_deepmrm[~vebios_deepmrm['peptide_id'].isin(odd_pept)].to_csv("Export/Vebios_ComplemEdge_deepmrm_20241108.csv", index=False)

# Write down peptides_counts
vebios_method.to_csv("export/peptide_counts.csv", index=True)

