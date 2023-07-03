# BCAC_Multi_ans_subtypes

- Goal  
  - To find out variants associate with breast cancer subtypes
  - To develop a PRS score for breast cancer subtypes
- Data  
It includes genotype and phenotype data of 95 studies from Breast Cancer Association Consortium ([BCAC](https://duckduckgo.com)). 
  - 40262 controls and 39784 cases from iCOGS array, and 68150 controls and 82314 cases from OncoArray
  - Genotype data were imputed using TOPMed (https://imputation.biodatacatalyst.nhlbi.nih.gov)
  - The data were stored on Biowulf at /data/BB_Bioinformatics/ProjectData/BCAC
  
- Code
  - check_data.R: Get the sample size tables and generate sample files to extract genotype data
  - read_gen.sh: Read and QC genotype data
  - form_genotype_block.sh: Form genotype data based on LD blocks
  - intrinsic_subtypes_genome.R: Run subype association analysis
  
