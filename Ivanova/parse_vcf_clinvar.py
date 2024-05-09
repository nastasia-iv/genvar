import csv
import gzip
import os
from typing import List

# csv.field_size_limit(sys.maxsize)

headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'CLNSIG', 'CLNVC', 'GENEINFO', 'MC',
           'Consequence', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'cDNA_position', 'CANONICAL']

# Column names for Clinvar VCF file
clinvar_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
clinvar_columns_dict = {col: idx for idx, col in enumerate(clinvar_columns)}

# Names for Clinvar VCF file (INFO column)
clinvar_info_names = [
    'AF_ESP', 'AF_EXAC', 'AF_TGP', 'ALLELEID', 'CLNDN', 'CLNDNINCL', 'CLNDISDB', 'CLNDISDBINCL',
    'CLNHGVS', 'CLNREVSTAT', 'CLNSIG', 'CLNSIGCONF', 'CLNSIGINCL', 'CLNVC', 'CLNVCSO', 'CLNVI',
    'DBVARID', 'GENEINFO', 'MC', 'ONCDN', 'ONCDNINCL', 'ONCDISDB', 'ONCDISDBINCL', 'ONC', 'ONCINCL',
    'ONCREVSTAT', 'ONCCONF', 'ORIGIN', 'RS', 'SCIDN', 'SCIDNINCL', 'SCIDISDB', 'SCIDISDBINCL',
    'SCIREVSTAT', 'SCI', 'SCIINCL'
]
clinvar_info_dict = {info: idx for idx, info in enumerate(clinvar_info_names)}

# Column names for VEP data
vep_names = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|'
             'HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|'
             'Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL').split('|')
vep_names_dict = {name: idx for idx, name in enumerate(vep_names)}


def parse_clinvar_vcf(vcf_file: str, output_dir: str = '', output_filename: str = '') -> None:
    """
    Parse a ClinVar VCF file and write the parsed data to a TSV file.

    Args:
        vcf_file (str): Path to the input ClinVar VCF file.
        output_dir (str, optional): Directory where the output TSV file will be saved. Defaults to current working directory.
        output_filename (str, optional): Name of the output TSV file. Defaults to the input VCF filename with '.tsv' extension.

    Returns:
        None
    """
    if output_dir == '':
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

    with gzip.open(vcf_file, 'rt') as input_file, open(output_file, 'a', newline='', encoding='utf-8') as output_file:
        vcf_reader = csv.reader(input_file, delimiter='\t')
        tsv_writer = csv.writer(output_file, delimiter='\t')

        for line in vcf_reader:
            parsed_data = parse_vcf_line(line)
            for data in parsed_data:
                tsv_writer.writerow(data)


def parse_vcf_line(line: List[str]) -> List[List[str]]:
    """
    Parse a single line of a VCF file.

    Args:
        line (List[str]): List representing a line from the VCF file.

    Returns:
        List[List[str]]: Parsed data for variant in the line.
    """
    parsed_data = []

    if not line[clinvar_columns_dict['CHROM']].startswith('#'):
        chrom = 'chr' + line[clinvar_columns_dict['CHROM']]
        position = line[clinvar_columns_dict['POS']]
        variation_id = line[clinvar_columns_dict['ID']]
        ref = line[clinvar_columns_dict['REF']]
        alt = line[clinvar_columns_dict['ALT']]

        info = line[clinvar_columns_dict['INFO']].split(';')
        for element in info:
            for category in clinvar_info_dict:
                if element.startswith(f'{category}='):
                    value = element.split('=')
                    clinvar_info_dict[category] = value[-1]

        vep_info = info[-1].split(',')

        for element in vep_info:
            variant_info = element.split('|')

            info_consequence = variant_info[vep_names_dict.get('Consequence')]
            info_symbol = variant_info[vep_names_dict.get('SYMBOL')]
            info_gene = variant_info[vep_names_dict.get('Gene')]
            info_feature_type = variant_info[vep_names_dict.get('Feature_type')]
            info_feature = variant_info[vep_names_dict.get('Feature')]
            info_biotype = variant_info[vep_names_dict.get('BIOTYPE')]
            info_cDNA = variant_info[vep_names_dict.get('cDNA_position')]
            info_canonical = variant_info[vep_names_dict.get('CANONICAL')]

            filtered_data = [
                chrom, position, variation_id, ref, alt,
                clinvar_info_dict['CLNSIG'], clinvar_info_dict['CLNVC'], clinvar_info_dict['GENEINFO'],
                clinvar_info_dict['MC'],
                info_consequence, info_symbol, info_gene, info_feature_type,
                info_feature, info_biotype,
                info_cDNA, info_canonical
            ]
            parsed_data.append(filtered_data)

    return parsed_data
