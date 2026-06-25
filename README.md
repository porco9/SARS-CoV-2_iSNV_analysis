# SARS-CoV-2_iSNV_analysis
Analysis of SARS-CoV-2 intra-host single nucleotide variants (iSNVs) under vaccine-induced immune pressure in a hamster model.
## Pipeline Overview

This repository contains a stepwise pipeline for identifying and analyzing intra-host single nucleotide variants (iSNVs) of SARS-CoV-2.

### Data processing workflow

1. Merge mapping coverage replicates  
   → `01_merge_mapping_coverage_replicates.py`

2. Apply read depth (≥100×) and minor allele frequency (5–50%) filters  
   → `02_filter_mapping_coverage_and_identify_iSNVs.py`

3. Generate iSNV summary table  
   → `03_build_iSNV_summary_table.py`

4. Remove iSNVs detected in the positive control group  
   → `04_remove_positive_control_iSNVs.py`

5. Map iSNVs to viral proteins and count variants  
   → `05_map_iSNVs_to_proteins_and_count.py`

6. Annotate synonymous and nonsynonymous mutations  
   → `06_annotate_synonymous_and_nonsynonymous_iSNVs.py`

7. Calculate callable sites per protein  
   → `07_calculate_protein_callable_sites.py`

### Statistical analysis

- iSNV rate and sharing analysis (Fig. 3)  
  → `Fig3_code.R`

- Negative binomial modeling and selection analysis (Fig. 4)  
  → `Fig4_code.R`
