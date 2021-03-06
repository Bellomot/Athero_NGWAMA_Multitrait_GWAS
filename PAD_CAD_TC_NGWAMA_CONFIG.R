########## Configuration File for PAD CAD TC scan

### Paths of summmary stats
Pathtrait1="/project/voight_GWAS/wbone/PAD_and_CAD_Bivariate_scan/CLEANEDforSHARE.MVP.EUR.PAD.results.anno_no_quotes.txt"

trait1_SNPcol = "ID"
trait1_EAcol = "EFFECT_ALLELE"
trait1_OAcol = "OTHER_ALLELE"
trait1_EAFcol = "EAF"
trait1_SEcol = "SE"
trait1_PVALcol = "PVAL"
trait1_BETAcol = "BETA"
trait1_BPcol = "POS"
trait1_CHRcol = "CHROM"
trait1_SAMPLESIZEcol = 243060

Pathtrait2="/project/voight_datasets/GWAS/21_CADMeta/CAD_META"

trait2_SNPcol = "oldID"
trait2_EAcol = "Allele1"
trait2_OAcol = "Allele2"
trait2_EAFcol = "Freq1"
trait2_SEcol = "StdErr"
trait2_PVALcol = "P-value"
trait2_BETAcol = "Effect"
trait2_BPcol = "BP"
trait2_CHRcol = "CHR"
trait2_SAMPLESIZEcol =547261

Pathtrait3="/project/voight_GWAS/wbone/NGWAMA/TC_remeta_data_with_MAF_Clean.txt"

trait3_SNPcol = "SNP"
trait3_EAcol = "A1"
trait3_OAcol = "A2"
trait3_EAFcol = "MAF"
trait3_SEcol = "SE"
trait3_PVALcol = "P"
trait3_BETAcol = "Beta"
trait3_BPcol = "BP"
trait3_CHRcol = "CHR"
trait3_SAMPLESIZEcol = 404128



### For each trait, add a trait lable to the traits data frame
traits <- as.data.frame(c("PAD","CAD","TC"))



### Create a list where there is one set of terms per trait of interest containing terms that will be used to search for existing GWAS catalog hits
strListTrait_total <- c(c("Peripheral artery", "Artery disease"), c("Heart Disease","Coronary","Artery disease"), c("Cholesterol","Total cho","LDL","HDL"))



### Data read in
# Summary stats
trait1 <- fread(Pathtrait1)
join1 <- trait1 %>% dplyr::select("chr.exposure"=trait1_CHRcol,"pos.exposure"=trait1_BPcol,"SNP"=trait1_SNPcol)
trait2 <- fread(Pathtrait2)
join2 <- trait2 %>% dplyr::select("chr.exposure"=trait2_CHRcol,"pos.exposure"=trait2_BPcol,"SNP"=trait2_SNPcol)
trait3 <- fread(Pathtrait3)
join3 <- trait3 %>% dplyr::select("chr.exposure"=trait3_CHRcol,"pos.exposure"=trait3_BPcol,"SNP"=trait3_SNPcol)

# 2 sample MR for NGWAMA
n_trait1 = read_exposure_data(Pathtrait1,sep="\t",snp_col=trait1_SNPcol,effect_allele_col=trait1_EAcol,other_allele_col=trait1_OAcol,eaf_col=trait1_EAFcol,se_col=trait1_SEcol,pval_col=trait1_PVALcol,beta_col=trait1_BETAcol)
n_trait1$samplesize.exposure <- trait1_SAMPLESIZEcol
n_trait1$samplesize.exposure <- as.numeric(n_trait1$samplesize.exposure)
n_trait2 = read_exposure_data(Pathtrait2,sep="\t",snp_col=trait2_SNPcol,effect_allele_col=trait2_EAcol,other_allele_col=trait2_OAcol,eaf_col=trait2_EAFcol,se_col=trait2_SEcol,pval_col=trait2_PVALcol,beta_col=trait2_BETAcol)
n_trait2$samplesize.exposure <- trait2_SAMPLESIZEcol
n_trait2$samplesize.exposure <- as.numeric(n_trait2$samplesize.exposure)
n_trait3 = read_exposure_data(Pathtrait3,sep="\t",snp_col=trait3_SNPcol,effect_allele_col=trait3_EAcol,other_allele_col=trait3_OAcol,eaf_col=trait3_EAFcol,se_col=trait3_SEcol,pval_col=trait3_PVALcol,beta_col=trait3_BETAcol)
n_trait3$samplesize.exposure <- trait3_SAMPLESIZEcol
n_trait3$samplesize.exposure <- as.numeric(n_trait3$samplesize.exposure)

# 2 sample MR for filtering
dat_trait1 = read_exposure_data(Pathtrait1,sep="\t",snp_col=trait1_SNPcol,effect_allele_col=trait1_EAcol,other_allele_col=trait1_OAcol,eaf_col=trait1_EAFcol,se_col=trait1_SEcol,pval_col=trait1_PVALcol,beta_col=trait1_BETAcol)
dat_trait2 = read_outcome_data(Pathtrait2,sep="\t",snp_col=trait2_SNPcol,effect_allele_col=trait2_EAcol,other_allele_col=trait2_OAcol,eaf_col=trait2_EAFcol,se_col=trait2_SEcol,pval_col=trait2_PVALcol,beta_col=trait2_BETAcol)
dat_trait3 = read_outcome_data(Pathtrait3,sep="\t",snp_col=trait3_SNPcol,effect_allele_col=trait3_EAcol,other_allele_col=trait3_OAcol,eaf_col=trait3_EAFcol,se_col=trait3_SEcol,pval_col=trait3_PVALcol,beta_col=trait3_BETAcol)



