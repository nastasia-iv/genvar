import csv
import os

import cyvcf2


def vcf_parsing(file_path: str) -> None:
    '''
    Parses the VCF file and extracts relevant data, then saves the processed data to a TSV file.

    Args:
    file_path (str): The path to the VCF file for parsing.

    Returns:
    str: A message indicating the creation of the TSV file.
    '''
    
    folder_path = 'processed_data'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f'Folder "{folder_path}" created.')
    else:
        print(f'Folder "{folder_path}" exists.')
    destination_file = f'processed_data/{file_path.split("/")
                                         [-1].split(".")[0]}.tsv'
    vcf = cyvcf2.VCF(file_path)
    print('VCF file downloaded')

    # Define the columns to extract
    info_fields_to_extract = ['AC', 'AC_afr', 'AC_amr', 'AC_nfe',
                              'AC_asj', 'AC_sas', 'AC_eas', 'AC_mid', 'AC_fin',
                              'AN', 'AN_afr', 'AN_amr', 'AN_nfe', 'AN_asj',
                              'AN_sas', 'AN_eas', 'AN_mid', 'AN_fin',
                              'AF', 'AF_afr', 'AF_amr', 'AF_nfe', 'AF_asj',
                              'AF_sas', 'AF_eas', 'AF_mid', 'AF_fin', 'vep']
    vep_field_mapping = {
            1: 'Consequence', 2: 'IMPACT', 3: 'SYMBOL', 4: 'Gene',
            5: 'Feature_Type', 6: 'Feature', 7: 'BIOTYPE', 8: 'EXON',
            9: 'INTRON', 17: 'ALLELE_NUM', 21: 'VARIANT_CLASS',
            24: 'CANONICAL', 44: 'LoF', 45: 'LoF_filter',
            46: 'LoF_flags', 47: 'LoF_info'}

    column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AC', 'AC_afr',
                    'AC_amr', 'AC_nfe', 'AC_asj', 'AC_sas', 'AC_eas',
                    'AC_mid', 'AC_fin', 'AN', 'AN_afr', 'AN_amr', 'AN_nfe',
                    'AN_asj', 'AN_sas', 'AN_eas', 'AN_mid', 'AN_fin', 'AF',
                    'AF_afr', 'AF_amr', 'AF_nfe', 'AF_asj', 'AF_sas', 'AF_eas',
                    'AF_mid', 'AF_fin', 'Consequence', 'IMPACT', 'SYMBOL',
                    'Gene', 'Feature_Type', 'Feature', 'BIOTYPE', 'EXON',
                    'INTRON', 'ALLELE_NUM', 'VARIANT_CLASS', 'CANONICAL',
                    'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']
    print('VCF parsing in progress...')

    # Initialize an empty list to store the extracted data
    data = []

    # Iterate over each variant in the VCF file
    for variant in vcf:
        if 'PASS' in variant.FILTERS:
            variant_data = [variant.CHROM, variant.POS,
                            variant.ID, variant.REF, variant.ALT[0]]
            info_data = [variant.INFO.get(field, '.')
                         for field in info_fields_to_extract]
            vep_annotation = variant.INFO.get('vep')

    # Handle multiple transcripts in vep if present
            if vep_annotation:
                vep_transcripts = vep_annotation.split(',')
                for transcript in vep_transcripts:
                    split_transcript = transcript.split('|')
                    vep_fields = []
                    for key in vep_field_mapping.keys():
                        try:
                            vep_fields.append(split_transcript[key])
                        except Exception:
                            vep_fields.append('.')
                    data.append(variant_data + info_data[:-1] + vep_fields)
            else:
                data.append(variant_data + info_data + ['.'])
    print('Data collected')
    with open(destination_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        # Write the header
        writer.writerow(column_names)
        # Write the data rows
        writer.writerows(data)
    return (f'{file_path.split("/")[-1].split(".")[0]}.tsv file created')
