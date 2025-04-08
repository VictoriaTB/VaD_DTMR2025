
######################################### 
### Taylor-Bateman et al., 2025 - Repurposing drugs for the prevention of vascular dementia: Evidence from drug target Mendelian randomization  

### Cis_DTMR.R requires 2 input files
### An Exposure GWAS 
### An Outcome GWAS

### Column should be labeled as follows
## ID - rsid of snp
## CHROM - Chromosome
## POS - Base pair location in build GRCh38_hg38
## Effect_allele - 
## Other_allele -
## Beta - Effect
## EAF - effect allele frequency if available
## SE - standard error
## P.value 

## The following information about your gene of interest is required 
Gene_Name <-       ##### Gene idenifer 
Gene_UB <-         ##### Gene location upper bound 
Gene_LB <-         ##### Gene location lower bound

### The folling information is needed on the outcome
Outcome_Name <-       ##### Outcome Name
  
library(TwoSampleMR)
library(dplyr)
library(plyr)
library(stringr)
library(ieugwasr)
library(MendelianRandomization)
### A access token is required
### Further information about the token system can be found here https://api.opengwas.io/api/ 

### The following outputs will be generated
### A table listing the IVs selected
### A table with the results of the MR analysis including IVW or wald ratio along with sensitivity analysis of MR-egger and weighted median so long as there are sufficient instruments 

#### Create output tables
Output_IV_Info <- data.frame(matrix(ncol = 16, nrow = 0))
x <- c("beta.exposure", "beta.outcome", "effect_allele.exposure",   
       "effect_allele.outcome", "exposure", "id.exposure", "other_allele.exposure", 
       "other_allele.outcome", "outcome", "pval.exposure", "pval.outcome",
       "SNP" , "SNP_index", "Gene", "Outcome", "Analysis" )
colnames(Output_IV_Info) <- x

Output_Results <- data.frame(matrix(ncol = 28, nrow = 0))
x <- c("Protein","Outcome", "Analysis",  
       "IVW_SNPs", "IVW_Est", "IVW_SE" , "IVW_lCI", "IVW_uCI", "IVW_Pvalue", "FStat", 
       "WM_SNPs", "WM_Est", "WM_SE" , "WM_lCI", "WM_uCI", "WM_Pvalue",
       "Egg_SNPs", "Egg_Est", "Egg_SE" , "Egg_lCI", "Egg_uCI", "Egg_Pvalue",
       "Egg_Int", "Egg_IntSE" , "Egg_IntlCI", "Egg_IntuCI", "Egg_IntPvalue"
)
colnames(Output_Results) <- x

####### Analysis

### subsetting the cis-acting gene region
### Set to include 500kb either side of CIS region

MaxV <- Gene_UB + 500000
LowV <- Gene_LB - 500000

Exposure_GWAS <-  subset(Exposure_GWAS, POS %in% c(paste(LowV):paste(MaxV) )) ### extract cis-acting gene region

Exposure_GWAS <- format_data(Exposure_GWAS, 
                            type = "exposure",
                            snp_col = "ID",
                            beta_col = "Beta",
                            se_col = "SE",
                            pval_col = "P.value",
                            chr_col = "CHROM",
                            pos_col = "POS",
                            effect_allele_col = "Effect_allele",
                            other_allele_col = "Other_allele",
                            eaf_col = "EAF")  

Exp_GWAS <-  subset(Exposure_GWAS, Exposure_GWAS$pval.exposure < 1e-04 ) ### subset to minimise data for clumping

### clumping, 10000kb region r2 < 0.001, genome wide significance
Clumping <- "5e-08"
Exp_SNPS <- clump_data(Exp_GWAS, clump_kb = 10000,
                        clump_r2 = 0.001,
                        clump_p1 = paste(Clumping),
                        clump_p2 = paste(Clumping),
                        pop= "EUR")

#### If no IVs found switches to reduced clumping threshold
if(nrow(Exp_SNPS) == 0) {
  Clumping <- "5e-05"
  Exp_SNPS <- clump_data(Exp_GWAS, clump_kb = 10000,
                          clump_r2 = 0.001,
                          clump_p1 = paste(Clumping),
                          clump_p2 = paste(Clumping),
                          pop= "EUR")}

Out_GWAS <- format_data(Outcome_GWAS, 
                          type = "outcome",
                          snp_col = "ID",
                          beta_col = "Beta",
                          se_col = "SE",
                          pval_col = "P.value",
                          chr_col = "CHROM",
                          pos_col = "POS",
                          effect_allele_col = "Effect_allele",
                          other_allele_col = "Other_allele",
                          eaf_col = "EAF") 


### Match Exp_SNPS
matched <- intersect( Out_GWAS$SNP, Exp_SNPS$SNP)
all <-  union(Out_GWAS$SNP, Exp_SNPS$SNP)
non.matched <- all[!all %in% matched]
Out_SNPS <- Out_GWAS[ which( ! Out_GWAS$SNP 
                             %in% non.matched) , ]

non.matched <- Exp_SNPS$SNP[!Exp_SNPS$SNP 
                            %in% matched]

Exp_SNPS <- Exp_SNPS[ which( Exp_SNPS$SNP 
                             %in% matched) , ]

Exp_IVs <- Exp_SNPS[order(Exp_SNPS$SNP),]
Out_IVs <- Out_SNPS[order(Out_SNPS$SNP),]

### Harmonise the two datasets
Harmonised_dat <- harmonise_data(exposure_dat = Exp_IVs, outcome_dat = Out_IVs)

#### Convert for use in MendelianRandomization Package
MR_Obj <- dat_to_MRInput(Harmonised_dat, get_correlations = FALSE, pop = "EUR")

MR1_pIV_pExp <- MR_Obj[[1]]

### Runs MR & generates results
new.row <- data.frame(Protein = Gene_Name, Outcome = Outcome_Name, Analysis = Clumping,
                      IVW_SNPs = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$SNPs) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$SNPs),
                      IVW_Est = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$Estimate) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$Estimate),
                      IVW_SE = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$StdError) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$StdError),
                      IVW_lCI = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$CILower) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$CILower),
                      IVW_uCI = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$CIUpper) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$CIUpper),
                      IVW_Pvalue = if (is.numeric(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$Pvalue) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw( MR1_pIV_pExp)$Pvalue),
                      FStat = if (is.numeric(MendelianRandomization::mr_ivw(MR1_pIV_pExp)$Fstat) == FALSE) print(NA) else print(MendelianRandomization::mr_ivw(MR1_pIV_pExp)$Fstat),
                      WM_SNPs = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$SNPs) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$SNPs),
                      WM_Est = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$Estimate) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$Estimate),
                      WM_SE = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$StdError) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$StdError),
                      WM_lCI = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$CILower) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$CILower),
                      WM_uCI = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$CIUpper) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$CIUpper),
                      WM_Pvalue = if (is.numeric(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$Pvalue) == FALSE) print(NA) else print(MendelianRandomization::mr_median( MR1_pIV_pExp, weighting = "weighted")$Pvalue),
                      Egg_SNPs = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$SNPs) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$SNPs),
                      Egg_Est = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Estimate) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Estimate),
                      Egg_SE = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$StdError.Est) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$StdError.Est),
                      Egg_lCI = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CILower.Est) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CILower.Est),
                      Egg_uCI = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CIUpper.Est) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CIUpper.Est),
                      Egg_Pvalue = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Pvalue.Est) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Pvalue.Est),
                      Egg_Int = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Intercept) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Intercept),
                      Egg_IntSE = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$StdError.Int) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$StdError.Int),
                      Egg_IntlCI = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CILower.Int) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CILower.Int),
                      Egg_IntuCI = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CIUpper.Int) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$CIUpper.Int),
                      Egg_IntPvalue = if (is.numeric(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Pvalue.Int) == FALSE) print(NA) else print(MendelianRandomization::mr_egger( MR1_pIV_pExp)$Pvalue.Int)
)

#### Adds to results table
Output_Results <- rbind(Output_Results, new.row)

### Saves harmonization table for reference
SNP_Data <- Harmonised_dat
SNP_Data <- SNP_Data %>% select(order(colnames(SNP_Data)))
SNP_Data$Gene <- paste(as.character(Gene_Name))
SNP_Data$Outcome <- paste(as.character(Outcome_Name))
SNP_Data$Analysis <- paste(as.character(Clumping))

SNP_Data<- SNP_Data[c("beta.exposure", "beta.outcome", "effect_allele.exposure",   
                      "effect_allele.outcome", "exposure", "id.exposure", "other_allele.exposure", 
                      "other_allele.outcome", "outcome", "pval.exposure", "pval.outcome",
                      "SNP" , "SNP_index", "Gene", "Outcome", "Analysis" )]

Output_IV_Info <- rbind(Output_IV_Info, SNP_Data)

#### This can be repeated for all outcomes and positive controls

