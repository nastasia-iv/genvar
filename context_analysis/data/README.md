## Data tables  
This directory contains tables obtained during sequence context analysis.  
  
* #### [all_nmd_escape_final_wo_in_clinvar.csv](intermediate_data/all_nmd_escape_final_wo_in_clinvar.csv)  
  CSV file with filtered stopgain variants (_pathogenic_ and _benign_) that escape nonsense-mediated decay (NMD). Variants marked as pathogenic in the Clinvar database were removed from the _benign_ dataset.  
    
  
* #### [all_nmd_undergo_final_wo_in_clinvar.csv](intermediate_data/all_nmd_undergo_final_wo_in_clinvar.csv)  
  CSV file with filtered stopgain variants (_pathogenic_ and _benign_) that undergo NMD. Variants marked as pathogenic in the Clinvar database were removed from the _benign_ dataset.  
     
         
* #### [clinvar_final_nmd_escape_df.csv](clinvar_final_nmd_escape_df.csv)  
  CSV file with pathogenic/likely pathogenic [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331) variants that escape NMD.   

    
* #### [clinvar_final_nmd_undergo_df.csv](clinvar_final_nmd_undergo_df.csv)  
  CSV file with pathogenic/likely pathogenic [Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) (v.20240331) variants falling under NMD.    

* #### [lof_loeuf_df.csv](intermediate_data/lof_loeuf_df.csv)  
  CSV file with stopgain variants located on canonical transcripts and not containing loss-of-function flags/filters according to the VEP annotation. Information obtained from exome [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) data (autosomes only). Variants are balanced by LOEUF score.  
    
* #### [lof_final_df.csv](lof_final_df.csv)  
  CSV file with stopgain variants located on canonical transcripts and not containing loss-of-function flags/filters according to the VEP annotation. Information obtained from exome [gnomad v4](https://gnomad.broadinstitute.org/downloads#v4) data (autosomes only). Variants are balanced by LOEUF and pext scores.  

    
* #### [nmd_escape_df.csv](nmd_escape_df.csv)  
  CSV file with filtered stopgain variants (gnomAD + Clinvar) that avoid NMD. Contains _no_ sequence context or codon information. Based on lof_final_df.csv.      

    
* #### [nmd_undergo_df.csv](nmd_undergo_df.csv)
  CSV file with filtered stopgain variants  (gnomAD + Clinvar) that undergo NMD. Contains _no_ sequence context or codon information. Based on lof_final_df.csv.    

