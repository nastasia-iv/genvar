import csv
import gzip
import os
from typing import List, Dict, Any

headers = [
    'Chr', 'Position', 'rsID', 'Ref', 'Alt', 'AC', 'Impact', 'Consequence',
    'Gene_symbol', 'Canonical_transcript', 'cDNA_position', 'LoF', 'LoF_flag', 'LoF_filter'
]

# Column names for VCF file
vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
vcf_columns_dict = {col: idx for idx, col in enumerate(vcf_columns)}

# Column names for VEP data
vep_names = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|'
            'cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|'
            'VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|'
            'UNIPROT_ISOFORM|SOURCE|||DOMAINS|miRNA|HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|'
            'MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info').split('|')

# Column names for population data
population_names = (
    'AC', 'AC_afr', 'AC_amr', 'AC_nfe', 'AC_asj', 'AC_sas', 'AC_eas', 'AC_mid', 'AC_fin',
    'AN', 'AN_afr', 'AN_amr', 'AN_nfe', 'AN_asj', 'AN_sas', 'AN_eas', 'AN_mid', 'AN_fin',
    'AF', 'AF_afr', 'AF_amr', 'AF_nfe', 'AF_asj', 'AF_sas', 'AF_eas', 'AF_mid', 'AF_fin'
)
population_dict = {pop: idx for idx, pop in enumerate(population_names)}


def parse_vcf(vcf_file: str, output_dir: str = '', output_filename: str = '') -> None:
    """
    Parse VCF file and write relevant data to a TSV file.

    Args:
        vcf_file (str): Path to the VCF file.
        output_dir (str): Directory to save the output TSV file. Default is current working directory.
        output_filename (str): Name of the output TSV file. Defaults to the input VCF filename with '.tsv' extension.

    Returns:
        None
    """
    if not output_dir:
        output_dir = os.getcwd()

    if output_filename == '':
        output_filename = os.path.splitext(os.path.basename(vcf_file))[0] + '.tsv'
    else:
        if not output_filename.endswith('.tsv'):
            output_filename += '.tsv'

    output_file = os.path.join(output_dir, output_filename)

    with open(output_file, 'w', newline='', encoding='utf-8') as table_file:
        writer = csv.writer(table_file, delimiter='\t')
        writer.writerow(headers)

    with gzip.open(vcf_file, 'rt') as input_file:
        vcf_reader = csv.reader(input_file, delimiter='\t')

        for line in vcf_reader:
            if line[vcf_columns_dict['CHROM']].startswith('chr') and line[vcf_columns_dict['FILTER']] == 'PASS':
                filtered_data = parse_line(line)
                write_to_output(output_file, filtered_data)


def parse_line(line: List[str]) -> List[str]:
    """
    Process a line from VCF file.

    Args:
        line (List[str]): List representing a line from the VCF file.

    Returns:
        List[str]: Processed data with selected information from the line.
    """
    chrom = line[vcf_columns_dict['CHROM']]
    position = line[vcf_columns_dict['POS']]
    rs_id = line[vcf_columns_dict['ID']]
    ref = line[vcf_columns_dict['REF']]
    alt = line[vcf_columns_dict['ALT']]
    population_dict.update(get_population_data(line))
    transcript_info = get_canonical_info(line)

    filtered_data = [
        chrom,
        position,
        rs_id,
        ref,
        alt,
        population_dict['AC'],
        ', '.join(transcript_info['IMPACT']),
        ', '.join(transcript_info['Consequence']),
        ', '.join(transcript_info['SYMBOL']),
        ', '.join(transcript_info['Feature']),
        ', '.join(transcript_info['cDNA_position']),
        ', '.join(transcript_info['LoF']),
        ', '.join(transcript_info['LoF_flags']),
        ', '.join(transcript_info['LoF_filter'])
    ]

    return filtered_data


def get_population_data(line: List[str]) -> Dict[str, str]:
    """
    Extract population data from VCF line.

    Args:
        line (List[str]): List representing a line from the VCF file.

    Returns:
        Dict[str, str]: Dictionary containing variant population data.
    """
    column_with_info = line[vcf_columns_dict['INFO']].split(';')
    frequency = {}
    for element in column_with_info:
        for pop in population_dict:
            if element.startswith(f'{pop}='):
                frequency[pop] = element.split('=')[-1]
    return frequency


def get_canonical_info(line: List[str]) -> Dict[str, List[str]]:
    """
    Extract canonical transcript information from VCF line.

    Args:
        line (List[str]): List representing a line from the VCF file.

    Returns:
        Dict[str, List[str]]: Dictionary containing information on canonical Ensemble transcript.
    """
    transcript_info = {field: [] for field in ['IMPACT', 'Consequence', 'SYMBOL', 'Feature', 'cDNA_position', 'LoF', 'LoF_flags', 'LoF_filter']}
    info = line[vcf_columns_dict['INFO']].split(';')
    vep_info = info[-1].split('|')
    feature_index = vep_names.index('Feature')
    canonical_index = vep_names.index('CANONICAL')
    symbol_index = vep_names.index('SYMBOL')
    step = 47

    for i in range(feature_index, len(vep_info), step):
        transcript = vep_info[i]
        canonical = vep_info[i + canonical_index - feature_index]
        symbol = vep_info[i + symbol_index - feature_index]

        if symbol and transcript.startswith('ENST') and canonical:
            for field in transcript_info:
                field_index = vep_names.index(field)
                info = vep_info[i + field_index - feature_index]
                if info:
                    if field == 'Feature' or field == 'CANONICAL':
                        transcript_info[field].append(transcript)
                    else:
                        transcript_info[field].append(info)
    return transcript_info


def write_to_output(output_file: str, data: List[str]) -> None:
    """
    Write data to the output TSV file.

    Args:
        output_file (str): Path to the output TSV file.
        data (List[str]): Data to be written to the file.

    Returns:
        None
    """
    with open(output_file, 'a', newline='', encoding='utf-8') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t')
        tsv_writer.writerow(data)
