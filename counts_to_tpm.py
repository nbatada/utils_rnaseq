import pandas as pd
import numpy as np
import warnings

def load_counts(counts_file_path):
    counts_df = pd.read_csv(counts_file_path, index_col=0)
    return counts_df

def load_gene_lengths(gene_lengths_file_path, gene_id_column):
    gene_lengths_df = pd.read_csv(gene_lengths_file_path, sep='\t')
    gene_lengths_subset = gene_lengths_df[[gene_id_column, 'Gene_Length']].copy()
    gene_lengths_subset = gene_lengths_subset.drop_duplicates(subset=[gene_id_column], keep='first')
    gene_lengths_series = gene_lengths_subset.set_index(gene_id_column)['Gene_Length']
    return gene_lengths_series

def filter_low_counts(counts_df, min_counts = 10, num_data = 3):
    filtered_df = counts_df[(counts_df >= min_counts).sum(axis=1) >= num_data]
    print('before:',counts_df.shape, 'after:',filtered_df.shape)
    return filtered_df

def align_data(counts_df, gene_lengths_series):
    aligned_counts_df = counts_df.reindex(gene_lengths_series.index).dropna(how='all')
    aligned_lengths = gene_lengths_series.reindex(aligned_counts_df.index)
    return aligned_counts_df, aligned_lengths

def compute_tpm(aligned_counts_df, aligned_lengths):
    aligned_lengths = aligned_lengths.replace(0, 1)
    aligned_counts_numeric = aligned_counts_df.apply(pd.to_numeric, errors='coerce')
    aligned_counts_numeric = aligned_counts_numeric.dropna(how='all')
    aligned_lengths = aligned_lengths.reindex(aligned_counts_numeric.index)
    aligned_lengths = aligned_lengths.replace(0, 1)
    rpk_df = aligned_counts_numeric.div(aligned_lengths / 1000, axis=0)
    scaling_factors = rpk_df.sum(axis=0, skipna=True)
    scaling_factors[scaling_factors == 0] = 1
    tpm_df = rpk_df.div(scaling_factors, axis=1) * 1e6
    return tpm_df

def plot_tpm_histograms(tpm_df, sampleinfo=None, facet_col='condition', sample_id_colname='sample_id'):
    import plotnine as p9
    tpm_df_melted = tpm_df.stack().reset_index()
    # Rename the column containing sample IDs to match sample_id_colname
    tpm_df_melted.columns = ['Gene', sample_id_colname, 'TPM']
    tpm_df_melted['log10_TPM_plus_1'] = np.log10(tpm_df_melted['TPM'] + 1)

    plot_data = tpm_df_melted
    if sampleinfo is not None and not sampleinfo.empty:
        # Merge using the specified sample_id_colname
        plot_data = pd.merge(tpm_df_melted, sampleinfo, on=sample_id_colname, how='left')

    plot = (
        p9.ggplot(plot_data, p9.aes(x='TPM', fill=facet_col)) # Use facet_col for fill
        + p9.geom_histogram(bins=100) # Removed fill from geom_histogram
        + p9.scale_x_log10()
        + p9.facet_wrap(f'~ {sample_id_colname} + {facet_col}', scales='fixed_y')
        + p9.labs(title='',
                  x='TPM',
                  y='Frequency')
        + p9.theme(figure_size=(10, 10))
    )

    return(plot)

# Main function
def main_counts_to_tpm(counts_df, gene_lengths_file_path='/Users/nbatada/lib/biolib/tpm/human_gene_lengths.tsv', gene_id_column='Gene_Symbol'):
    gene_lengths_series = load_gene_lengths(gene_lengths_file_path, gene_id_column)
    aligned_counts_df, aligned_lengths = align_data(counts_df, gene_lengths_series)
    tpm_df = compute_tpm(aligned_counts_df, aligned_lengths)
    return tpm_df

counts_file='counts.csv'
gene_lengths_file='human_gene_lengths.tsv'
counts_df = load_counts(counts_file)
tpm_df=main_counts_to_tpm(counts_df)


sampleinfo_df=pd.read_csv('sampleinfo_detailed.csv')
print(sampleinfo_df.columns)
plot_tpm_histograms(tpm_df, sampleinfo_df, facet_col='tissue_detailed',sample_id_colname='sample_id')

    
