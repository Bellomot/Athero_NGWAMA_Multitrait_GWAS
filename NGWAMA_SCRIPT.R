########## NGWAMA CODE





#Load the following on an interactive session before using this script:
#module load JAGS
#module load OpenBUGS
#module load plink/1.90Beta4.5
#module load bedtools2/2.25.0
#module load ldsc/1.0.1
#module load R/3.6.3
#source activate ldsc 

#Download the following data files:
#  plink version 1.9 (https://www.cog-genomics.org/plink2)
#  1000G eur reference files (http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz)

library(TwoSampleMR)
library(vroom)
library(data.table)
library(furrr)
library(glue)
library(ggpubr)
library(magrittr)
library(httr)
library(jsonlite)
library(fs)
library(broom)
library(tidyverse)
library(phenoscanner)
library(gt)
library(janitor)
library(ieugwasr)
library(tidygenomics)
library(ggQQunif)
library(ggupset)
library(miscF)

print((.packages()))

source("/project/voight_GWAS/wbone/NGWAMA/N_WEIGHTED_GWAMA_SOURCE_CODE.R")
source("NGWAMA_CONFIG.R")





#################################
##### Formatting for NGWAMA #####
#################################
#Create ts to use in future naming for written out files 
ts=paste0(traits[1,1],sep="")
for (i in 2:nrow(traits)){
  assign("ts", paste0(ts,traits[i,1],sep=""))
}

# Select columns and gzip to munge folder
for(i in 1:nrow(traits)) { 
assign(paste("n_trait",i,sep=""), dplyr::left_join(eval(parse(text=paste("n_trait",i,sep=""))),eval(parse(text=paste("join",i,sep=""))), by="SNP"))
eval(parse(text=paste("n_trait",i,sep=""))) %>% vroom_write(paste0("./munge/trait",i,".txt.gz"))
}





##################
##### NGWAMA #####
##################

### Munge sumstats
### must change directory here ###
# write munge script
dir_ls("./munge") %>%
  enframe() %>%
  mutate(name = str_replace(path_file(value), "\\.txt.gz$", "")) %>%
  mutate(munge_sumstats = glue("munge_sumstats.py --sumstats {value} --snp SNP --N-col samplesize.exposure --signed-sumstats beta.exposure,0 --a1 effect_allele.exposure --a2 other_allele.exposure --p pval.exposure --frq eaf.exposure --merge-alleles /project/voight_GWAS/wbone/NGWAMA/ldsc/ldsc-master/w_hm3.snplist --out ./munge_final/{name}")) %>%
  dplyr::select(munge_sumstats) %>%
  write_delim("./munge_traits.sh", col_names = FALSE, delim = "\t")
# run munge script
system(paste0("chmod 777 ./munge_traits.sh")) 
system(paste0("./munge_traits.sh")) 

### h2 analysis
# Write h2 script
dir_ls("./munge_final", regexp = ".gz$") %>%
  enframe() %>%
  mutate(name = str_replace(path_file(value), ".sumstats.gz", "")) %>%
  mutate(h2 = glue("ldsc.py --h2 {value} --ref-ld-chr /project/voight_GWAS/wbone/NGWAMA/ldsc/ldsc-master/eur_w_ld_chr/ --w-ld-chr /project/voight_GWAS/wbone/NGWAMA/ldsc/ldsc-master/eur_w_ld_chr/ --out ./{name}_h2")) %>%
  dplyr::select(h2) %>%
  write_delim("./h2_traits.sh", col_names = FALSE, delim = "\t")
# run h2 script
system(paste0("chmod 777 ./h2_traits.sh"))
system(paste0("./h2_traits.sh"))

### Cross trait ldsc
# Write cross-trait ldsc script
dir_ls("./munge_final", regexp = ".gz$") %>%
  crossing(p1 = ., p2 = .)  %>% 
  mutate(key = paste0(pmin(p1, p2), pmax(p1, p2), sep = "")) %>% 
  arrange(p1) %>%
  filter(p1 != p2) %>%
  distinct(key, .keep_all = TRUE) %>%
  dplyr::select(-key) %>%
  mutate(
    p1_name = str_replace(path_file(p1), "\\.sumstats.gz$", ""),
    p2_name = str_replace(path_file(p2), "\\.sumstats.gz$", "")
  ) %>%
  mutate(ldsc_correlation = glue("ldsc.py --rg {p1},{p2} --ref-ld-chr /project/voight_GWAS/wbone/NGWAMA/ldsc/ldsc-master/eur_w_ld_chr/ --w-ld-chr /project/voight_GWAS/wbone/NGWAMA/ldsc/ldsc-master/eur_w_ld_chr/ --out ./{p1_name}-{p2_name}_res")) %>%
  dplyr::select(ldsc_correlation) %>%
  write_delim("./ldsc_traits.sh", col_names = FALSE, delim = "\t")
# run cross trait script
system(paste0("chmod 777 ./ldsc_traits.sh"))
system(paste0("./ldsc_traits.sh"))

### Formatting
ldsc_res <- dir_ls("./", regexp = "\\_res.log$") %>%
  enframe(name = NULL) %>%
  mutate(ldsc = map(value, read_table2, skip = 60, n_max = 1)) %>%
  unnest(ldsc) %>%
  ungroup %>%
  dplyr::select(-value) %>%
  mutate(p1 = path_file(p1),
         p2 = path_file(p2)) %>%
  mutate(p1 = str_replace_all(p1, ".sumstats.gz", ""),
         p2 = str_replace_all(p2, ".sumstats.gz", ""))
###!! always specify the delim in vroom function as " "!!###
h2_res <- dir_ls("./", regexp = "\\_h2.log$") %>%
  enframe(name = NULL) %>%
  mutate(ldsc = map(value, vroom, skip = 25, n_max = 1, delim = " ", col_names = FALSE)) %>%
  unnest(ldsc) %>%
  ungroup %>%
  dplyr::select(value, X5) %>%
  mutate(h2 = str_replace(X5, " *\\(.*", "")) %>%
  mutate(value = str_replace(path_file(value), "_h2.log", "")) %>%
  dplyr::select(-X5) %>%
  type_convert()

cti_matrix <- ldsc_res %>%
  bind_rows(ldsc_res %>% rename(p1 = p2, p2 = p1)) %>%
  reshape2::acast(p1~p2, value.var = "gcov_int", fill = 1)

### Run NGWAMA script
# Create an empty list of however many triats there are
gwama_dat <- vector(mode = "list", length = nrow(traits))
# Use a repeat loop to insert items into the list 
x <- 1
repeat {
  gwama_dat[[x]] <- eval(parse(text=paste("n_trait",x,sep=""))) %>%
    mutate(Z = beta.exposure/se.exposure) %>%
    dplyr::select(SNPID = SNP, CHR = chr.exposure, BP = pos.exposure, EA = effect_allele.exposure, OA = other_allele.exposure, EAF = eaf.exposure, N = samplesize.exposure, Z, P = pval.exposure) %>%
    as.data.frame()

  x = x + 1
  if(x == nrow(traits)+1){
    break
  }
}
# NGWAMA script 
multivariate_GWAMA(x = gwama_dat, cov_Z = cti_matrix, h2 = h2_res$h2, out="./", name = "NGWAMA", output_gz = TRUE, check_columns = FALSE)

### Format the output file for filtering
ngwama <-vroom("./NGWAMA.N_weighted_GWAMA.results.txt.gz")
colnames(ngwama) <- c("CPTID","SNP","CHR","BP","effect_allele","other_allele","eaf.trait","maf.trait","samplesize.trait1","samplesize.trait2","direction","beta.trait","se.trait","z_trait","P")
ngwama <- as.data.table(ngwama)
write.table(ngwama, quote = FALSE,row.names=FALSE,file=paste0("AllSNPs_",ts,"_multivarresults_chisquared_pvalue.txt"),sep="\t")





###########################################
##### Create functions for processing #####
###########################################

# Create a function to search for strings of words related to each trait for pulling out significant SNPs
noMatch <-function(currTraits, strList_total){
  for (tStr in strList_total){
    if (any(grep(tStr,currTraits,ignore.case=TRUE)))
    {
      return(FALSE)
    }
  }
  return(TRUE)
} 

# Create a function to extract all GWAS hits from GWAS ctalog within 500kb of significant SNP without considering LD
findNearGWASCatalogHits<-function(r){
  Chrom = trimws(r["CHR"])
  Coord=r["BP"]
  res = c("","","","")
  col = paste0("chr",Chrom)
  for (i in 1:dim(GWASlist[[col]])[1]){
    G = GWASlist[[col]][i,]
    if (abs(as.numeric(Coord)-as.numeric(G["coord"]))<500000){
      res[1]= paste(res[1],G["trait"])
      res[2] = paste(res[2],G["PMID"])
      res[3] = paste(res[3],G["rs"])
      res[4] = paste(res[4],abs(as.numeric(Coord)-as.numeric(G["coord"])))
    }}
  return (res)
}

# Create a function to see if any SNP proxies (???) are in the GWAS catalog
f = function(SNPProxies){
  leadSNP = SNPProxies["SNP"]
  SNPProxyList = unlist(strsplit(as.character(SNPProxies["TAGS"]),'|',fixed=TRUE))
  matches = which(GWAS[,"rs"] %in% SNPProxyList)
  if(length(matches)==0){
    return(c(leadSNP,"NONE","NONE","NONE"))
  }
  return(c(leadSNP,paste(GWAS[matches,"rs"],collapse = ","),paste(GWAS[matches,"trait"],collapse = ","),paste(GWAS[matches,"PMID"],collapse=",")))
}

# Create a function to Find the most significant SNP in a block of SNPs using the single trait GWAS p-values
findHighestSNP = function(proxySnps,GWASFile1,pvalcolName){
  proxyLines = GWASFile1[which(GWASFile1$SNP %in% unlist(proxySnps)),]
  maxLine = proxyLines[which.min(unlist(proxyLines[pvalcolName])),]
  return(maxLine)
}

# Create a function to extract the lowest p value SNP in each haploblock from each trait GWAS
getProxies = function(SNPLine){
  SNPProxyList = append(SNPLine["SNP"],strsplit(as.character(SNPLine["TAGS"]),'|',fixed=TRUE))
  return(SNPProxyList)
}

# Exit functin if no novel SNPs are found
exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}  





#####################################################################################
##### Harmonize the data just so that I can use this data set for a merge later #####
#####################################################################################

### Harmonize data frames and change column names
x <- 2
repeat {
  assign(paste("har",x,sep=""), harmonise_data(exposure_dat=dat_trait1, outcome_dat= eval(parse(text=paste("dat_trait",x,sep="")))))
  
  har_colnames <- colnames(eval(parse(text=(paste("har",x,sep="")))))
  har_colnames <- gsub("exposure","trait1", har_colnames)
  har_colnames <- gsub("outcome",paste("trait",x,sep=""), har_colnames)
  
  setnames(eval(parse(text=(paste("har",x,sep="")))), c(har_colnames))
  
  x = x + 1
  if(x == nrow(traits)+1){
    break
  }
}

### Merge Harmonized data frames base on all exposure (trait1) columns
#To repeat loop merge, you must first merge two data frames to have a common data frame
if(nrow(traits)>2){
  har_total <- merge(har2, har3, by = c("SNP","effect_allele.trait1","other_allele.trait1","beta.trait1","eaf.trait1","se.trait1","pval.trait1","trait1","pval_origin.trait1"))
}

x <- 4
if(nrow(traits)>3){
  repeat {
    har_total <- merge(har_total, eval(parse(text=paste("har",x,sep=""))), by = c("SNP","effect_allele.trait1","other_allele.trait1","beta.trait1","eaf.trait1","se.trait1","pval.trait1","trait1","pval_origin.trait1"))
    
    x = x + 1
    if(x == nrow(traits)+1){
      break
    }
  }
}

if(nrow(traits)==2){
  har_total <- har2
}

### Take the highest p value per SNP out of all traits and z score p values (otherwise, will have underflow AKA the p values may be too small to be reported by the computer)
for(i in 1:nrow(traits)){
  har_total[,paste("pval.trait",i,sep="")] <- pmax(har_total[,paste("pval.trait",i,sep="")],1e-16)
}

# Calculate a z score for the p value of each trait for N variate scan slater 
for(i in 1:nrow(traits)){
  har_total[,(ncol(har_total)+1)] <- sqrt(qchisq(1 - har_total[,paste("pval.trait",i,sep="")],1))*sign(har_total[,paste("beta.trait",i,sep="")])
  colnames(har_total)[ncol(har_total)] <- paste("z_trait",i,sep="")
}

# Select columns of interest and write out data 
write.table(har_total[,grepl("SNP|effect_allele.trait|other_allele.trait|eaf.trait|beta.trait|pval.trait|se.trait|z_trait",colnames(har_total))],paste0("AllSNPs_",ts,"_multivarpvalues_zscores.txt"),quote = FALSE,row.names=FALSE)




#####################################################
##### Format output file from Multivariate Scan #####
#####################################################

### Perform LD-clumping of Multivariate Results to  get clumps/loci of significant SNPs 
system(paste0("plink --noweb --bfile /project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed  --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink --clump-p1 .000001 --clump-r2 0.2  --clump-kb 1000 --clump AllSNPs_",ts,"_multivarresults_chisquared_pvalue.txt --out allChrMergedClumped"))
cData = read.table("allChrMergedClumped.clumped",header=TRUE,stringsAsFactors = FALSE)



### Find the gene nearest to the singificant SNP and create a bed format from a clumped file 
system("awk ' {OFS=\"\t\"} (NR>1 && NF>0) {print \"chr\"$1,$4,$4+1,$3}' allChrMergedClumped.clumped | sort -k1,1 -k2,2n>  allChrMergedClumped.bed")  #Puts top SNP into bedfile format
system(paste0("closestBed -a allChrMergedClumped.bed -b /project/voight_selscan/ksiewert/BetaGenomeScan/genes_biomart_sorted.bed -d > ",ts,"_NearestGene.bed"))

# Format nearest gene output to merge with original SNP file
nearGene = read.table(paste0(ts,"_NearestGene.bed"),header=FALSE,stringsAsFactors = FALSE)
colnames(nearGene) = c("chr","coord1","coord2","SNP","chr2","c1","c2","gene","dist")
nearGene = aggregate(nearGene[,c("gene","dist")],list(nearGene$SNP),toString)
colnames(nearGene)[1] = "SNP"

# Merge nereast gene file with significant SNP data set 
d = merge(cData,har_total[,grepl("SNP|effect_allele.trait|other_allele.trait|eaf.trait|beta.trait|se.trait|pval.trait|z_trait", colnames(har_total))],by="SNP")
d = merge(d,nearGene[c("SNP","gene","dist")],by="SNP")




### Filter for loci that are nominally significant
o = order(d$P,decreasing=FALSE)
d = d[o,]
d = d[which(d$P<1e-7),] #contains all clustered SNP with p<1e-7 whether or not novel
write.table(d, quote = FALSE,row.names=FALSE,file=paste0("AllSigSNPs_",ts,"_r2clumped.txt"),sep="\t") #NOTE: May not be 1MB from each other, but will be r2<.2





#############################################################
###### Select only Novel Significant SNPs from SNP list #####
#############################################################

### Output all hits from the GWAS catalog within 500kb of nominally significant SNPs for any trait 
GWAS = read.table("/project/voight_datasets/GWAS/15_GWASCatalog/Jan_2019_version/GWAShg37.bed",header=FALSE,sep="\t",colClasses = "character",fill=TRUE,quote = "",stringsAsFactors = FALSE)
# Rename GWAS catalog columns
colnames(GWAS)=c("chr","coord","coord+1","PMID","rs","trait","P")
# Filter GWAS for significant SNPS only
GWAS = GWAS[which(as.numeric(GWAS$P)<=5e-8),] 
GWASlist = split(GWAS,GWAS$chr)
# Apply findNearGWASCatalogHits function created earlier to extract GWAS hits within 500kb of singificant SNP
GWASoverlap = t(apply(d,1,findNearGWASCatalogHits)) 
colnames(GWASoverlap) = c("nearby_trait","nearby_PMID","nearby_rs","nearby_coord")

### Filter out Significant SNPs that are already annotated by GWAS results 
# Bind singificant SNP dataset and GWAS overlap output 
GWASannot = cbind(d,GWASoverlap)
write.table(GWASannot, quote = FALSE,row.names=FALSE,file=paste0("AllSigSNPs_",ts,"_r2clumped_GWAScat_annotated.txt"),sep="\t")

# Use noMatch function to search the GWAS catalog hits for strings of words related to each trait and return all SNPs not annotated by GWAS
g = unlist(lapply(GWASoverlap[,1],noMatch,strListTrait_total)) 
# Select SNPs without any GWAS catalog hits
novel = GWASannot[which(g),]
# Select SNPs with p-values from each trait GWAS that are nominally, but not genome-wide significant (if were genome wide significant in each trait GWAS, then would have already been reported)  
pval_total <- colnames(novel)[grepl("pval.trait", colnames(novel))]
novelSNPs <- novel %>% dplyr::filter_at(vars(pval_total), all_vars(.>5e-8 & .<0.005))
# Write out for plink function 
write.table(novelSNPs$SNP,row.names=FALSE,col.names=FALSE,quote=FALSE,file="Multivar_novelSNPrsNums.txt",sep="\t")

### If there are no novel SNPs, exit script
if(nrow(novelSNPs)==0){
  print("No novel loci detected with this scan")
  print("The following traits were run:")
  print(ts)
  print("The GWAS catalog used:")
  print("/project/voight_datasets/GWAS/15_GWASCatalog/Jan_2019_version/GWAShg37.bed")
  exit()
}


###############################################################################
###### AND Select traits within r^2=.2 (reasonable LD) of novel GWAS hits #####
###############################################################################

# Compare Nvariate scan novel SNP to find other SNPS with R^2 0.2 in realtion to our novel SNPs
system(paste0("plink --bfile /project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink --show-tags Multivar_novelSNPrsNums.txt --out NovelSNPTags  --tag-kb 1000 --list-all --tag-r2 0.2"))
SNPs = read.table("NovelSNPTags.tags.list", header=TRUE,stringsAsFactors = FALSE)
overlappingTraits = t(apply(SNPs,1,f))
colnames(overlappingTraits)=c("SNP","Linked_SNPs","Linked_Trait","Linked_Link")
write.table(overlappingTraits,row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste0("SigSNPs_",ts,"_no_GWAScat_hits_500kb_r2_GWAScat_annotated.txt"),sep="\t")

# need to make sure none of the proxies are relevant GWAS hits-similar process to above but with r^2 instead of distance
# Use noMatch function to search the GWAS catalog hits for strings of words related to each trait and filter out 
notKnowns = lapply(overlappingTraits[,"Linked_Trait"],noMatch,strListTrait_total)
novelSNPs = novelSNPs[which(unlist(notKnowns)),]
novelSNPs = merge(novelSNPs,overlappingTraits,by="SNP")
save(novelSNPs,file="novelSNPs")
write.table(novelSNPs,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("NovelSNPs_",ts,"_no_GWAScat_hits_500kb_r2_GWAScat_filtered.txt"),sep="\t")

### Use getProxies function to select the lowest p value SNP in each haploblock and add haplogroup information (findHighestSNP)
for(i in 1:nrow(traits)) {
  dat_colnames <- colnames(eval(parse(text=paste("dat_trait",i,sep=""))))
  dat_colnames <- gsub("exposure","trait1", dat_colnames)
  dat_colnames <- gsub("outcome",paste("trait",i,sep=""), dat_colnames)
  setnames(eval(parse(text=paste("dat_trait",i,sep=""))), c(dat_colnames))
  
  assign(paste("top_trait",i,sep=""), t(sapply(apply(SNPs,1,getProxies),findHighestSNP,eval(parse(text=paste("dat_trait",i,sep=""))), paste("pval.trait",i,sep=""))))
  
  top_trait <- eval(parse(text=paste("top_trait",i,sep="")))
  colnames(top_trait)[colnames(top_trait)=='SNP'] <- paste("SNP_top_trait",i,sep="")
  assign(paste("top_trait",i,sep=""), top_trait)
}

### To make sure output of top_trait is a data frame, we first make a df, make all columns charater vectors, and then make df
for(i in 1:nrow(traits)) {
  top_trait <- eval(parse(text=paste("top_trait",i,sep="")))
  top_trait <- as.data.frame(top_trait)
  top_trait <- apply(top_trait,2,as.character)
  top_trait <- as.data.frame(top_trait)
  if(ncol(top_trait)==1) {
    top_trait <- t(top_trait)
    top_trait <- as.data.frame(top_trait)
  }
  assign(paste("top_trait",i,sep=""), top_trait)
}

### Bind back SNP columns for lowest p value SNPs in each haploblock and select columns of interest 
x <- 1

repeat{
  assign(paste("top_trait",x,sep=""), cbind(eval(parse(text=paste("top_trait",x,sep=""))), SNP=SNPs$SNP)) 
  
  assign(paste("top_trait",x,sep=""), eval(parse(text=paste("top_trait",x,sep="")))[,grepl("SNP|SNP_top_trait|beta.trait|se.trait|pval.trait",colnames(eval(parse(text=paste("top_trait",x,sep="")))))])
  
  x = x+1
  if(x == nrow(traits)+1) {
    break
  }
}

### Merge all outputs from lowest p value SNP in each haploblock back onto original data 
# Make novelSNPs a data frame so we can grepl
# Make novelSNPs a data frame so we can grepl
novelSNPs <- as.data.frame(novelSNPs[,grepl("SNP|CHR|BP|P|effect_allele.trait|other_allele.trait|eaf.trait|beta.trait|se.trait|pval.trait|z_trait|gene|dist|Direction|Linked_SNPs|Linked_Trait|Linked_Link|nearby_trait|nearby_PMID|nearby_rs", colnames(novelSNPs))])
# Merge top_trait1 and then start repeat loop
results = merge(novelSNPs,top_trait1,by="SNP")
write.table(results,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("NovelSNPs_",ts,"_results.txt"),sep="\t")
write.table(top_trait1,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("NovelSNPs_",ts,"_top_trait1.txt"),sep="\t")

# Repeat loop
x <- 2
repeat{
  results <- merge(results,eval(parse(text=paste("top_trait",x,sep=""))),by="SNP")
  
  x = x+1
  if(x == nrow(traits)+1) {
    break
  }
}

### Subsitute actual trait names for traiti
results_colnames <- colnames(results)
for(i in 1:nrow(traits)){
  results_colnames <- gsub(paste("trait",i,sep=""),traits[i,1], results_colnames)
}
colnames(results) <- results_colnames

# Write out results
fwrite(results,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("NovelSNPs_",ts,"_fully_annotated.txt"),sep="\t")




##################################################
##### Regional association plots using plink #####
##################################################

# Function to extract table of LD corrleation information (r)

ld_extract <- function(variants, bfile, plink_bin) {
  # Make textfile
  shell <- ifelse(Sys.info()["sysname"] == "Mac", "cmd", "sh")
  fn <- tempfile()
  
  data.frame(variants) %>%
    vroom::vroom_write(fn, col_names = FALSE)
  
  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --extract ", shQuote(fn, type = shell),
    " --r inter-chr ",
    " --out ", shQuote(fn, type = shell)
  )
  system(fun2)
  
  res <- data.table::fread(paste0(fn, ".ld")) %>%
    dplyr::select(SNP_A, SNP_B, R)
  
  system(paste0("rm ", fn, ".ld"))
  
  return(bind_rows(res, res %>% rename(SNP_A = SNP_B, SNP_B = SNP_A)))
}

# Function to create regional association plots using plink for clumping
gg_regional_association_plink <- function(df, lead_snps = NULL, rsid = rsid, chromosome = chromosome, position = position, p_value = p_value, p_value_threshold = 1, clump_kb = 1000, clump_r2 = 1, plot_distance = 500000, bfile, plink_bin, plot_title = NULL, plot_subtitle = NULL, n_row = 2) {
  df <- df %>%
    dplyr::select(rsid = {{ rsid }}, chromosome = {{ chromosome }}, position = {{ position }}, p_value = {{ p_value }}) %>%
    mutate_if(is.factor, as.character)
  
  if (!is.null(lead_snps)) {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(rsid %in% lead_snps)
  } else {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(pval < p_value_threshold) %>%
      ld_clump(bfile = bfile, plink_bin = plink_bin, clump_kb = 500, clump_r2 = 0.2)
  }
  
  locus_snps <- df %>%
    filter(rsid %in% indep_snps$rsid) %>%
    dplyr::select(chromosome, position, lead_rsid = rsid) %>%
    pmap_dfr(function(chromosome_filter = first, position_filter = second, lead_rsid = third) {
      df %>%
        filter(chromosome == chromosome_filter & between(position, position_filter - plot_distance / 2, position_filter + plot_distance / 2)) %>%
        mutate(lead_position = position_filter) %>%
        mutate(lead_rsid = lead_rsid)
    }) %>%
    mutate(lead_marker = glue::glue("{chromosome}:{lead_position}")) %>%
    group_by(lead_marker) %>%
    mutate(lead_p_value = min(p_value)) %>%
    ungroup() %>%
    mutate(label = paste0(lead_marker, "\n", lead_rsid)) %>%
    mutate(label = fct_reorder(label, lead_p_value))
  
  
  locus_snps_ld <- locus_snps %>%
    group_nest(lead_rsid, keep = TRUE) %>%
    ungroup() %>%
    mutate(r2_matrix = map(data, function(x) {
      ld_extract(x$rsid, bfile = bfile, plink_bin = plink_bin) %>%
        rename(lead_rsid = SNP_A, rsid = SNP_B) %>%
        filter(lead_rsid %in% indep_snps$rsid) %>%
        mutate(r2 = abs(R)^2) %>%
        dplyr::select(-R)
    })) %>%
    dplyr::select(r2_matrix) %>%
    unnest(r2_matrix) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    bind_rows(tibble(
      lead_rsid = unique(locus_snps$lead_rsid),
      rsid = unique(locus_snps$lead_rsid),
      r2 = rep(1, length(unique(locus_snps$lead_rsid)))
    )) %>%
    mutate(color_code = case_when(
      rsid == lead_rsid ~ "purple",
      TRUE ~ color_code
    )) %>%
    mutate(color_code = fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "blue", "blue4")) %>%
    mutate(lead = lead_rsid == rsid) 
  
  plot <- locus_snps %>%
    left_join(locus_snps_ld) %>%
    mutate(r2 = ifelse(is.na(r2) == TRUE, 0.1, r2)) %>%
    mutate(lead = ifelse(is.na(lead) == TRUE, FALSE, lead)) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    mutate(color_code = ifelse(lead == TRUE, "purple", color_code)) %>%
    ggplot(aes(position / 1000000, -log10(p_value), fill = color_code, size = lead, alpha = lead, shape = lead)) +
    geom_point() +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
    scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = c("Lead SNP", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"), breaks = c("purple", "red", "orange", "darkgreen", "blue", "blue4")) +
    scale_size_manual(values = c(3, 8), guide = FALSE) +
    scale_shape_manual(values = c(21, 23), guide = FALSE) +
    scale_alpha_manual(values = c(0.8, 1), guide = FALSE) +
    scale_x_continuous(n.breaks = 3) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 8))) +
    facet_wrap(~label, scales = "free", nrow = n_row) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Position (Mb)",
      y = expression(-log[10]("p-value"))
    ) +
    theme_light(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title.align = 0.5,
      legend.key = element_rect(size = 6, fill = "white", colour = NA)
    )
  
  return(plot)
}



### Creating and saving plots
###!!Need a print() statement in order to save file as pdf!!###
# Using trait data loaded in from config file to select columns to pipe into gg_regional_association_plink function 

x <- 1

repeat{
  
  for (i in 1:nrow(results)) {
    
    rm(graph_df)
    rm(plot)
    
    graph_df <- eval(parse(text=paste("trait",x,sep="")))
    
    graph_df <- graph_df %>% dplyr::select(rsid = eval(parse(text=paste0("trait",x,"_SNPcol"))), chromosome = eval(parse(text=paste0("trait",x,"_CHRcol"))), position = eval(parse(text=paste0("trait",x,"_BPcol"))), p_value = eval(parse(text=paste0("trait",x,"_PVALcol"))))
    
    head(graph_df)
    
    plot <-  gg_regional_association_plink(graph_df,p_value_threshold = 5e-8, lead_snps = results[i,1] , bfile = "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref", plink_bin = "plink", plot_title = paste0("NGWAMA trait",x," ",print(results[i,1])," Regional Association Plot"), plot_subtitle = expression("Bivariate Variants p < 1x10^-7"))
    
    pdf(file = paste0("./TRAIT",x,"_",print(results[i,1]),"_RA_PLOTS.pdf"), paper = 'USr', width = 15, height = 20)
    print(plot)
    dev.off()
  }
  
  x = x + 1
  if (x == nrow(traits)+1) {
    break
  }
}

