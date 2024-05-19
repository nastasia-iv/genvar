from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def get_context(df: pd.DataFrame, transcript_fasta: dict, left_len: int, right_len: int) -> str:
    """
    Add sequences context to the dataframe.

    Args:
        df (pd.DataFrame): The dataframe containing transcript and position information.
        transcript_fasta (dict): A dictionary containing transcript sequences.
        left_len (int): The length of the left flanking region.
        right_len (int): The length of the right flanking region.

    Returns:
        str: A message indicating that contexts have been added to the dataframe.
    """
    contexts = []

    for index, row in df.iterrows():
        transcript_id = row['Canonical_transcript']
        position_of_interest = row['cDNA_position']

        if transcript_id in transcript_fasta:
            transcript_length = len(transcript_fasta[transcript_id])

            if left_len < position_of_interest < transcript_length - right_len:
                sequence_of_interest = str(transcript_fasta[transcript_id][position_of_interest - left_len: position_of_interest + right_len])
            else:
                sequence_of_interest = None
        else:
            sequence_of_interest = None
        
        contexts.append(sequence_of_interest)

    df['Context'] = contexts
    message = 'Contexts have been added to the dataframe!'

    return message


def check_ref(row: pd.Series, transcript_fasta: dict) -> str:
    """
    Checks if the reference matches the target letter in the context. 
    At the same time, it determines whether the variant is on the + or - strand.

    Args:
        row: pd.Series: A row from a DataFrame.
        transcript_fasta: dict: A dictionary containing transcript sequences.

    Returns:
        str: A string indicating whether the reference base matches the variant in the genomic context.
            Returns '+' strand if they match, '-' strand if the variant is complementary to the reference, and
            'Not_defined' if the transcript ID is not found or other conditions are not met.
    """
    transcript_id = row['Canonical_transcript']
    cDNA_position = row['cDNA_position']
    variant = str(row['Context'][12])
    ref = row['REF']

    complement_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    if transcript_id in transcript_fasta:
        fasta_sequence = str(transcript_fasta[transcript_id][int(cDNA_position) - 1])

        if variant == fasta_sequence and variant == ref:
            return '+'
        elif variant == fasta_sequence == fasta_sequence and variant == complement_bases.get(ref):
            return '-'
        else:
            return 'Not_defined'
    else:
        return 'Not_defined'


def get_codon_info(row: pd.Series) -> Union[int, str]:
    """
    Change the reference option to an alternative one in the context of the sequence. 
    Then determine at which codon position the variant is located and which stop codon results.

    Args:
        row: pd.Series: A row from a DataFrame.

    Returns:
        Union[Tuple[int, str], str]: A tuple containing the codon position and the codon itself, or a string indicating the absence of a stop codon or strand information.
    """
    context = row['Context']
    strand = row['Strand']
    alt = row['ALT']
    complement_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    if strand == '+':
        updated_context = context[:12] + alt + context[13:]
    elif strand == '-':
        complement_alt = complement_bases.get(alt)
        updated_context = context[:12] + complement_alt + context[13:]
    else:
        return 'No_strand', 'No_strand'

    codon_position_1 = updated_context[10:13]
    codon_position_2 = updated_context[11:14]
    codon_position_3 = updated_context[12:15]

    if codon_position_1 in ['TAA', 'TAG', 'TGA']:
        return 3, codon_position_1
    elif codon_position_2 in ['TAA', 'TAG', 'TGA']:
        return 2, codon_position_2
    elif codon_position_3 in ['TAA', 'TAG', 'TGA']:
        return 1, codon_position_3
    else:
        return 'No_stop', 'No_stop'



def filter_and_convert_to_list(column: pd.Series) -> List[str]:
    """
    Filter and convert a dataframe column to a list.

    Args:
        column (pd.Series): The dataframe column to filter and convert.

    Returns:
        List[str]: A list containing filtered values from the dataframe column.
    """
    filtered_list = column.tolist()
    filtered_list = list(filter(None, filtered_list))
    return filtered_list


def calculate_chi2_p_values(context_ben: (List[str]), context_pat: (List[str])) -> Tuple[List[float], List[float]]:
    """
    Calculate chi-squared values and p-values (with FDR) for two sets of context sequences.

    Args:
        context_ben (List[str]): List of context sequences for benign samples.
        context_pat (List[str]): List of context sequences for pathogenic samples.

    Returns:
        Tuple[List[float], List[float]]: Chi-squared values and p-values.
    """
    sequences_array_ben = np.array(context_ben)
    sequences_array_pat = np.array(context_pat)

    freq_array_ben = np.zeros((len(sequences_array_ben[0]), 4))
    freq_array_pat = np.zeros((len(sequences_array_pat[0]), 4))

    for seq in sequences_array_ben:
        for j, nt in enumerate(seq):
            if nt == 'A':
                freq_array_ben[j][0] += 1
            elif nt == 'C':
                freq_array_ben[j][1] += 1
            elif nt == 'G':
                freq_array_ben[j][2] += 1
            elif nt == 'T':
                freq_array_ben[j][3] += 1

    for seq in sequences_array_pat:
        for j, nt in enumerate(seq):
            if nt == 'A':
                freq_array_pat[j][0] += 1
            elif nt == 'C':
                freq_array_pat[j][1] += 1
            elif nt == 'G':
                freq_array_pat[j][2] += 1
            elif nt == 'T':
                freq_array_pat[j][3] += 1

    chi2_values = []
    p_values = []

    for i in range(freq_array_pat.shape[0]):
        contingency_table = np.array([freq_array_pat[i], freq_array_ben[i]])
        try:
            chi2, p_value, _, _ = chi2_contingency(contingency_table)
        except ValueError as e:
            # pseudo p-value in case of error
            chi2_values.append(np.nan)
            p_values.append(1.0)
        else:
            chi2_values.append(chi2)
            p_values.append(p_value)

    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

    return chi2_values,  p_values_corrected
