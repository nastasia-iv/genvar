from typing import List

import numpy as np
import pandas as pd


def info_collecting(input_files: List[str]) -> pd.DataFrame:
    '''
    Collects and combines information from
    multiple input TSV files into a single DataFrame.

    Args:
    input_files (List[str]): List of input file names.

    Returns:
    pd.DataFrame: Combined DataFrame
        containing the information from all input files.
    '''

    combined_df = pd.read_csv(input_files[0], sep='\t', low_memory=False)
    for file in input_files[1:]:
        df = pd.read_csv(file, sep='\t', low_memory=False, header=None)
        df.columns = combined_df.columns
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df


def info_filtering(gene_data: pd.DataFrame,
                   constraint_file: str,
                   expression_file: str) -> pd.DataFrame:
    '''
    Filters and analyzes genomic data for specific information related to
    gene expression and constraints.

    Args:
    gene_data (pd.DataFrame): DataFrame containing genomic data.
    constraint_file (str): File path for constraint data.
    expression_file (str): File path for expression data.

    Returns:
    pd.DataFrame: Processed DataFrame
        containing the filtered and analyzed genomic information.
    '''

    population_ac = ['AC_afr', 'AC_amr', 'AC_nfe', 'AC_asj',
                     'AC_sas', 'AC_eas', 'AC_mid', 'AC_fin']
    transcript_list = []
    sum_ac = []
    sum_population_ac = [[], [], [], [], [], [], [], []]
    gene_names = []
    gene_id = []
    alt_sum = []
    loeuf_values = []
    exon_numbers = []
    expression = []

    # Reading data
    constraint_transcript = pd.read_csv(constraint_file, sep='\t')
    expression_transcript = pd.read_csv(expression_file, sep='\t')

    # Filtration dataframe
    gene_data = gene_data[(gene_data['Feature_Type'] == 'Transcript') &
                          (gene_data['BIOTYPE'] == 'protein_coding')]
    values_to_filter = ['stop_gained', 'frameshift_variant',
                        'splice_donor_variant', 'splice_acceptor_variant']
    gene_data = gene_data[
        gene_data['Consequence'].isin(values_to_filter)
        ][gene_data['Feature'].str.contains('ENST')]

    columns_to_keep = ['gene', 'gene_id', 'transcript', 'canonical',
                       'lof.oe_ci.upper', 'num_coding_exons']
    constraint_transcript_loeuf = constraint_transcript[columns_to_keep]

    # Collecting specific data
    sum_ac_per_transcript = gene_data.groupby('Feature')['AC']
    for key, group in sum_ac_per_transcript:
        transcript_list.append(key)
        unique_values = group.sum()
        sum_ac.append(unique_values)

    for idx, el in enumerate(population_ac):
        sum_population = []
        sum_AC_per_transcript = gene_data.groupby('Feature')[el]
        for key, group in sum_AC_per_transcript:
            unique_values = group.sum()
            sum_population.append(unique_values)
            sum_population_ac[idx].append(sum(sum_population))
            sum_population = []

    gene_name_per_transcript = gene_data.groupby('Feature')['SYMBOL']
    for key, group in gene_name_per_transcript:
        unique_values = group.unique()
        gene_names.extend(unique_values)

    gene_id_per_transcript = gene_data.groupby('Feature')['Gene']
    for key, group in gene_id_per_transcript:
        unique_values = group.unique()
        gene_id.extend(unique_values)

    alt_per_transcript = gene_data.groupby('Feature')['ALT']
    for key, group in alt_per_transcript:
        unique_values = len(group.sum())
        alt_sum.append(unique_values)

    max_ac_indices = gene_data.groupby('Feature')['AC'].idxmax()
    result = gene_data.loc[max_ac_indices, ['Feature', 'Consequence', 'AC']]
    freq_cons = list(result['Consequence'])
    freq = list(result['AC'])

    # New dataframe grouping
    transcripts_df_ac = pd.DataFrame({
        'Transcript_ID': transcript_list,
        'AC': sum_ac,
        'AC_afr': sum_population_ac[0],
        'AC_amr': sum_population_ac[1],
        'AC_nfe': sum_population_ac[2],
        'AC_asj': sum_population_ac[3],
        'AC_sas': sum_population_ac[4],
        'AC_eas': sum_population_ac[5],
        'AC_mid': sum_population_ac[6],
        'AC_fin': sum_population_ac[7],
        'Gene_name': gene_names,
        'Gene_id': gene_id,
        'Variant': alt_sum,
        'Max_AC_in_transcript': freq,
        'Consequence_of_max_AC': freq_cons
        })

    # AC/N metric
    transcripts_df_ac['AC/Variant'] = transcripts_df_ac['AC'] /\
        transcripts_df_ac['Variant']

    # LOEUF metric and exon number
    transcripts = transcripts_df_ac['Transcript_ID'].unique()
    tr = []
    for transcript in transcripts:
        # LOEUF metric and exon number
        loeuf = list(constraint_transcript[constraint_transcript['transcript'].str.contains(transcript)]['lof.oe_ci.upper'])
        exon = list(constraint_transcript[constraint_transcript['transcript'].str.contains(transcript)]['num_coding_exons'])
        # Expression level
        expression_level = list(expression_transcript[expression_transcript['ID_transcript'].str.contains(transcript)]['Max_median_expression'])
        if loeuf != []:
            loeuf_values.append(loeuf[0])
        else:
            loeuf_values.append(np.nan)
        if exon != []:
            exon_numbers.append(exon[0])
        else:
            exon_numbers.append(np.nan)
        if expression_level != []:
            expression.append(expression_level[0])
        else:
            expression.append(np.nan)

    transcripts_df_ac['LOEUF_transcript'] = loeuf_values
    transcripts_df_ac['Exon_number'] = exon_numbers
    transcripts_df_ac['Max_median_expression'] = expression

    return transcripts_df_ac
