# Athero_NGWAMA_Multitrait_GWAS
# Multi-trait GWAS code for Bellomo et al "Multi-trait GWAS of atherosclerosis detects novel pleiotropic loci"
All of the code, configuration files, and dependency files to perform an NGWAMA multi-trait GWAS and eQTL colocalization using GTEx v8 are in this repository. 

Instructions on running the NGWAMA multi-trait GWAS pipeline: 
  - Create an analysis directory in which all outputs from the script will be saved. In that analysis directory, have the following items:
    - Save the desired NGWAMA multi-trait GWAS script (NGWAMA_SCRIPT.R) and edit this script to point to the location of the  "N_WEIGHTED_GWAMA_SOURCE_CODE.R". 
    - Save the desired configuration file (NGWAMA_CONFIG.R) corresponding to the GWAS files you would like to perform a multi-trait GWAS on. 
    - Create two empty directories named "munge" and "munge_final" which will be used by the pipeline.  
  - Lastly, execute the NGWAMA multi-trait GWAS in the analysis directory (eg Rscript ./NGWAMA_SCRIPT.R). 

Code file:
  - NGWAMA_SCRIPT.R => This script will take in any n number of GRCh37 GWAS summary statistic files with the following information for each SNP: rs number corresponding to the SNP, effect allele, alternate allele, effect allele frequency, standard error, P value, beta, position in GRCh37, and chromosome. These files will be used to perform a multi-trait NGWAMA scan to detect novel loci associated with all traits. 
  - Known errors in this file include early exit of the script due to absence of novel SNPs. To ensure this is the case, check the log files for the printout "No novel loci detected with this scan". 

Individual configuration files: the configuration files here were used to perform 15 different scans of combinations of atherosclerotic traits. NOTE: if you want to replicate our analysis, rename the configuration file (NGWAMA_CONFIG.R), change the directories to the locations of the GRCh37 GWAS summary statistic files, and follow the instructions above for running the NGWAMA multi-trait GWAS pipeline above. 
  - PAD_BMI_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and body mass index bivariate GWAS. 
  - PAD_CAD_BMI_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and body mass index trivariate GWAS. 
  - PAD_CAD_HDL_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and high density lipoprotein trivariate GWAS. 
  - PAD_CAD_LDL_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and low density lipoprotein trivariate GWAS. 
  - PAD_CAD_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and coronary artery disease bivariate GWAS.
  - PAD_CAD_SMK_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and smoking trivariate GWAS. 
  - PAD_CAD_T2D_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and type 2 diabetes trivariate GWAS. 
  - PAD_CAD_TC_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and total cholesterol trivariate GWAS. 
  - PAD_CAD_TG_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease, coronary artery disease, and triglyceride trivariate GWAS. 
  - PAD_HDL_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and high density lipoprotein bivariate GWAS. 
  - PAD_LDL_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and low density lipoprotein bivariate GWAS. 
  - PAD_SMK_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and smoking bivariate GWAS. 
  - PAD_T2D_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and type 2 diabetes bivariate GWAS. 
  - PAD_TC_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and total cholesterol bivariate GWAS. 
  - PAD_TG_NGWAMA_CONFIG.R => The configuration file for peripheral artery disease and triglyceride bivariate GWAS. 

Dependency files:
  - N_WEIGHTED_GWAMA_SOURCE_CODE.R => source code for NGWAMA script downloaded from (https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/N_weighted_GWAMA.function.1_2_6.R?raw=TRUE).
  - data_maf0.01_rs_ref => 1000 Genome European reference files downloaded from (http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz). 
  - genes_biomart_sorted.bed => gene positions bedfile downloaded from biomart (http://www.biomart.org).  
  - EUR.final => A file containing the 503 population IDs associated with 1KG Phase III vcf, which was self-generated. 
  - EUR.final.plink => A plink file than contains data for only the individuals from EUR.final. 
  - GWAShg37.bed => BED file from the January 2019 GWAS catalog flat file converted into a BED file. 


  
Instructions on running the eQTL colocalization pipeline with files in the (EQTL_COLOC_CONFIG) directory:
  - EQTL Colocalizer pipeline information can be found on <brian github link>. In brief, create an analysis directory in which individual folders for each SNP output will be saved. In that analysis directory, please save the following items:
    - Save the desired eQTL colocalization script (eqtl_colocalizer.R). 
    - Save the desired configuration file (eQTL_config_template.R) corresponding to the GWAS file you would like to use for colocalizatin with eQTL form GTEx v8.
    - Save the shell script (bash_eQTL_COLOC_loop.sh) that contains a loop to create a folder and run the eqtl_colocalizer.R script for each SNP you wish to colocalize with eQTL data.  
- Lastly, execute the shell script in the analysis directory (eg ./bash_eQTL_COLOC_loop.sh). 
  
Code file:
  - EQTL Colocalizer full code can be found on <brian github link>. 
  - eqtl_colocalizer.R => This script will take a list of SNPs (eg CAD_15_COMBO_ALL_LEADSNPS.txt or PAD_15_COMBO_ALL_LEADSNPS.txt) with the following information: chromosome, start position in hg37, stop position in hg37, and rs number corresponding to the SNP. These SNPs will be extracted from one summary statistic file defined in the configuration file with the following information: rs number corresponding to the SNP, effect allele, alternate allele, effect allele frequency, standard error, P value, beta, position in hg37, chromosome, and number of participants. The information from the sumary statistics will be colocalized with eQTL data from GTEx v8. 
  - Known errors in this file include if the filename field for an eGene-Tissue pair is NA. To ensure this is the case, check the log files for the printout "This tissue is not available in the all pairs files currently". 

Individual configuration file: the configuration files here were used to run 2 sets of colocalization based on either peripheral artery disease or coronary artery disease sumary statistics. NOTE: if you want to replicate our analysis, rename the R script based on the desired trait (eg PAD_eqtl_colocalizer.R), change the directory to the locations of the hg37 GWAS summary statistic files, and follow the instructions above for running the eQTL colocalization pipeline above. 
  - eQTL_config_template.R => The configuration file for peripheral artery disease summary statistics and eQTL GTEx v8 colocalization. 

Dependency files:
  - data_maf0.01_rs_ref => 1000 Genome european reference files downloaded from (http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz). 
  - GTEx_v8_Tissue_Summary_with_filenames.csv => A file of GTEx single trait eQLTs from the GTEx online GUI. More information can be found on <brian github link>. 
  - GTEx_Analysis_v8_eQTL_tabix => Tabix files for all SNPs downloaded from GTEx. More information can be found on <brian github link>. 
