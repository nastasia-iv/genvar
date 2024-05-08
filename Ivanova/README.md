## Brief content information
This directory contains code to perform the task of analyzing the effects of sequence context on variant frequencies. A list of required packages and libraries is presented in [requirements.txt](requirements.txt).

* ### [all_code.py](all_code.py)
  A set of functions for processing data: obtaining the sequence context, calculating chi-square and p-values when comparing contexts, and drawing plots of the dependence of p-values on the position of the variant.

  
* ### [data_processing.ipynb](data_processing.ipynb)  
  Jupiter notebook for processing data from of [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) exomes and [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331). It contains a brief data analysis, primary variant filtering, and dataframes create options. The resulting dataframes and images are in the [data](data) and [images](images) folders, respectively.

  
* ### [get_context_nmd_escape.ipynb](get_context_nmd_escape.ipynb)  
  Jupiter notebook with sequence context analysis for nonsense-mediated mRNA decay (NMD) escaped region variants. Graphs with analysis results are located in the [images](images). 

  
* ### [get_context_nmd_undergo.ipynb](get_context_nmd_undergo.ipynb)  
  Jupiter notebook with sequence context analysis for variants falling under NMD. Graphs with analysis results are located in the [images](images).  

  
* ### [parse_vcf_canonical.py](parse_vcf_canonical.py)  
  Function for obtaining information about [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) variants located on canonical Ensemble transcripts.

  
* ### [parse_vcf_clinvar.py](parse_vcf_clinvar.py)  
  Function for parsing pre-processed with VEP [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331) file.
