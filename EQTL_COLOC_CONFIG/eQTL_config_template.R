#ID for user's trait of interest. (Can be any string)
trait = "PAD"
#path to the input files
traitFilePath = "/project/voight_GWAS/wbone/PAD_and_CAD_Bivariate_scan/CLEANEDforSHARE.MVP.EUR.PAD.results.anno_no_quotes.txt"
#column IDs from trait file
trait_A1col = "EFFECT_ALLELE"
trait_A2col = "OTHER_ALLELE"
trait_SNPcol = "ID"
trait_CHRcol = "CHROM"
trait_BPcol = "POS"
trait_Pcol = "PVAL"
trait_Ncol = "N"
trait_MAFcol = "EAF"
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType = "cc"
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp = 0.1288  #look this up
#locus information for running coloc. Currently this assumes these genomic positions to be from hg19
chrom = CHROMOSOME
colocStart = STARTBP
colocStop = STOPBP
lead_SNP = "SNPNUMBER" 

