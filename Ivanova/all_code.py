import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency
from Bio.Seq import Seq


def get_context(df, transcript_fasta, left_len, right_len):
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


def filter_and_convert_to_list(column):
    filtered_list = column.tolist()
    filtered_list = list(filter(None, filtered_list))
    return filtered_list


def calculate_chi2_p_values(context_ben, context_pat):
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
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
        chi2_values.append(chi2)
        p_values.append(p_value)

    return chi2_values, p_values

def plot_p_values(positions_list, p_values_list):
    plt.figure(figsize=(10, 3))
    #colors = ['#C8A2C8', '#B0E0E6', '#FFB6C1']
    colors = ['lightcoral', 'lightblue', 'plum']

    for positions, p_values, color in zip(positions_list, p_values_list, colors):
        plt.plot(positions, np.log10(p_values), marker='o', linestyle='-', color=color)

    plt.title('Logarithm of P-Values Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Log10(p-value)')
    plt.grid(False)
    plt.legend(['1st', '2nd', '3rd'], title='Position in codon')

    all_positions = [pos for sublist in positions_list for pos in sublist]
    plt.xticks(all_positions, [str(pos) if pos == 0 else f"{pos}" if pos < 0 else f"+{pos}" for pos in all_positions])

    plt.axvspan(-0.5, 0.5, color='yellow', alpha=0.3)

    plt.show()


def translate_context(df_column):
    aa_column = []
    for context in df_column:
        if context is None:
            aa_column.append(None)
        else:
            sequence = Seq(context)
            aa_sequence = sequence.translate()
            aa_column.append(str(aa_sequence))
    return aa_column
