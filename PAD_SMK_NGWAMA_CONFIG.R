########## Configuration File for PAD SMK scan, note SMK is space seperated 

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

Pathtrait2="/project/voight_GWAS/wbone/NGWAMA/SMK_GWAS_UKBB/2019.10.02_Lifetime_Smoking_GWAS_Data_Sheet_1.txt"

trait2_SNPcol = "SNP"
trait2_EAcol = "EFFECT_ALLELE"
trait2_OAcol = "OTHER_ALLELE"
trait2_EAFcol = "EAF"
trait2_SEcol = "SE"
trait2_PVALcol = "P"
trait2_BETAcol = "BETA"
trait2_BPcol = "BP"
trait2_CHRcol = "CHR"
trait2_SAMPLESIZEcol = 462690



### For each trait, add a trait lable to the traits data frame
traits <- as.data.frame(c("PAD","SMK"))



### Create a list where there is one set of terms per trait of interest containing terms that will be used to search for existing GWAS catalog hits 
strListTrait_total <- c(c("Peripheral artery", "Artery disease"), c("Smoke","Smoking"))



### Data read in
# Summary stats
trait1 <- fread(Pathtrait1)
join1 <- trait1 %>% dplyr::select("chr.exposure"=trait1_CHRcol,"pos.exposure"=trait1_BPcol,"SNP"=trait1_SNPcol)
trait2 <- fread(Pathtrait2)
join2 <- trait2 %>% dplyr::select("chr.exposure"=trait2_CHRcol,"pos.exposure"=trait2_BPcol,"SNP"=trait2_SNPcol)

# 2 sample MR for NGWAMA
n_trait1 = read_exposure_data(Pathtrait1,sep="\t",snp_col=trait1_SNPcol,effect_allele_col=trait1_EAcol,other_allele_col=trait1_OAcol,eaf_col=trait1_EAFcol,se_col=trait1_SEcol,pval_col=trait1_PVALcol,beta_col=trait1_BETAcol)
n_trait1$samplesize.exposure <- trait1_SAMPLESIZEcol
n_trait1$samplesize.exposure <- as.numeric(n_trait1$samplesize.exposure)
n_trait2 = read_exposure_data(Pathtrait2,sep=" ",snp_col=trait2_SNPcol,effect_allele_col=trait2_EAcol,other_allele_col=trait2_OAcol,eaf_col=trait2_EAFcol,se_col=trait2_SEcol,pval_col=trait2_PVALcol,beta_col=trait2_BETAcol)
n_trait2$samplesize.exposure <- trait2_SAMPLESIZEcol
n_trait2$samplesize.exposure <- as.numeric(n_trait2$samplesize.exposure)

# 2 sample MR for filtering
dat_trait1 = read_exposure_data(Pathtrait1,sep="\t",snp_col=trait1_SNPcol,effect_allele_col=trait1_EAcol,other_allele_col=trait1_OAcol,eaf_col=trait1_EAFcol,se_col=trait1_SEcol,pval_col=trait1_PVALcol,beta_col=trait1_BETAcol)
dat_trait2 = read_outcome_data(Pathtrait2,sep=" ",snp_col=trait2_SNPcol,effect_allele_col=trait2_EAcol,other_allele_col=trait2_OAcol,eaf_col=trait2_EAFcol,se_col=trait2_SEcol,pval_col=trait2_PVALcol,beta_col=trait2_BETAcol)




