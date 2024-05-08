## Brief content information
This directory contains code to perform the task of analyzing the effects of sequence context on variant frequencies.

* ### [all_code.py](all_code.py)
  A set of functions for processing data: obtaining the sequence context, calculating chi-square and p-values when comparing contexts, drawing a plots of the dependence of p-values on the position of the variant.

  
* ### data_processing.ipynb
  Jupiter notebook for processing data from of gnomad v4 and Clinvar (v.20240331). Contains a brief data analysis, primary variant filtering and dataframes create options.

  
* ### get_context_nmd_escape.ipynb
  Jupiter notebook with sequence context analysis for variants in the nonsense-mediated mRNA decay (NMD) escape region. 

  
* ### get_context_nmd_undergo.ipynb
  Jupiter notebook with sequence context analysis for variants falling under NMD.

  
* ### parse_vcf_canonical.py  
  Function for obtaining information about gnomad v4 variants located on canonical Ensemble transcripts.

  
* ### parse_vcf_clinvar.py  
  Function for parsing pre-processed with VEP Clinvar (v.20240331) file.
