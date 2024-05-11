## 💫 Sequence context analysis
This directory contains code for analyzing the effect of sequence context on variant frequency. The list of required packages and libraries is presented in the [environment.yaml](environment.yaml). Please note: scripts for identifying NMD(-) regions must be obtained separately from [here](https://github.com/rebeccaito/nmd-escape).

* #### [analysis_functions.py](analysis_functions.py)
  A set of functions for processing data: obtaining the sequence context, calculating chi-square and p-values when comparing contexts, and drawing plots of the dependence of p-values on the position of the variant.

  
* #### [data_processing.ipynb](data_processing.ipynb)  
  Jupyter notebook for processing data from of [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) exomes and [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331). It contains a brief data analysis, primary variant filtering, and dataframes create options. The resulting dataframes and images are in the [data](data) and [images](images) folders, respectively.

  
* #### [get_context_nmd_escape.ipynb](get_context_nmd_escape.ipynb)  
  Jupyter notebook with sequence context analysis for nonsense-mediated mRNA decay (NMD) escaped region variants. Graphs with analysis results are located in the [images](images). 

  
* #### [get_context_nmd_undergo.ipynb](get_context_nmd_undergo.ipynb)  
  Jupyter notebook with sequence context analysis for variants falling under NMD. Graphs with analysis results are located in the [images](images).  

  
* #### [parse_vcf_canonical.py](parse_vcf_canonical.py)  
  Function for obtaining information about [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) variants located on canonical Ensemble transcripts.
  To run parser, import function `parse_vcf` as shown below, specifying the path to the compressed (`.bgz`) vcf file. If necessary, you can specify the output folder and file name. More details can be found in the function docstring.
  ```python
  from parse_vcf_canonical import parse_vcf  
  parse_vcf("path/to/file.bgz")  
  ```

  
* #### [parse_vcf_clinvar.py](parse_vcf_clinvar.py)  
  Function for parsing pre-processed with VEP [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331) file.
    To run Clinvar parser, import function `parse_clinvar_vcf` as shown below, specifying the path to the compressed (`.gz`) vcf file. If necessary, you can specify the output folder and file name. More details can be found in the function docstring.
  ```python
  from parse_vcf_clinvar import parse_clinvar_vcf  
  parse_clinvar_vcf('path/to/file.vcf.gz')
  ```
