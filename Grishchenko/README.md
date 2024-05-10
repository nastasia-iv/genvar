## ✨ Transcript Conservation Assessment

### Overview
This repository contains scripts and notebooks for the analysis of transcript conservation in human genes. 
The analysis combines information about genetic variations in transcripts with expression data and
classic conservativity metrics to gain insights into transcript conservation across genes. A list of required packages and libraries is presented in [environment.yaml](environment.yaml).

#### Code Directory
[**`vcf_parser.py`**](vcf_parser.py) - Script for parsing VCF files to collect necessary data.

> Usage:
> ```python
> from vcf_parser import vcf_parsing
> vcf_parsing('../data/NAME_OF_YOUR_FILE.vcf')
> ```

[**`gnomad_vcf_parser.ipynb`**](gnomad_vcf_parser.ipynb) - Jupyter notebook demonstrating the usage of `vcf_parsing` from vcf_parser.py.

[**`all_chr_genvar_analysis.ipynb`**](all_chr_genvar_analysis.ipynb) - Jupyter notebook illustrating gene analysis from all chromosomes,
focusing on chromosome 13.
The analysis integrates information on transcript variants, expression levels obtained from GTEx portal, and
conservativity metrics from gnomAD database.

[**`all_genes_analysis.ipynb`**](all_genes_analysis.ipynb) - Jupyter notebook providing insights into all human genes
and outlining the gene selection process for further investigation.

[**`expression_data_collecting.ipynb`**](expression_data_collecting.ipynb) - Script for collecting transcript expression data from GTEx Transcript TPMs data.

[**`gene_examples_analysis.ipynb`**](gene_examples_analysis.ipynb) - Analysis of transcript examples from selected genes for subsequent phases of the study.

#### Data directory

[**example_chr22.vcf**](example_chr22.vcf) - VCF file containing genetic variation data for chromosome 22.

[**max_tissue_median_expr.tsv**](max_tissue_median_expr.tsv) - TSV file with maximum median tissue expression data.

[**all_genes_analysis.tsv**](all_genes_analysis.tsv) - TSV file with detailed analysis of all human genes.

[**CFLAR_analysis.tsv**](CFLAR_analysis.tsv) - TSV file with analysis results for the gene CFLAR.

[**RPS10_analysis.tsv**](RPS10_analysis.tsv) - TSV file with analysis results for the gene RPS10.

[**HLADRA_analysis.tsv**](HLADRA_analysis.tsv) - TSV file with analysis results for the gene HLA-DRA.

[**TPT1_analysis.tsv**](TPT1_analysis.tsv) - TSV file with analysis results for the gene TPT1.

[Plots directory](https://github.com/nastasia-iv/BI_project_2024/tree/main/Grishchenko/data/plots) contains plots generated from the Jupyter notebooks included in the repository.

For additional details, refer to the specific files and directories within this repository.
