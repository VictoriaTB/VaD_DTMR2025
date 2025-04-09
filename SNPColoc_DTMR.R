######################################### 
### Taylor-Bateman et al., 2025 - Repurposing drugs for the prevention of vascular dementia: Evidence from drug target Mendelian randomization  

### SNPColoc_DTMR.R requires 2 input files
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
## N - Total sample size column 

## The following information about your gene of interest is required 
Gene_Name <-       ##### Gene idenifer 
Gene_UB <-         ##### Gene location upper bound 
Gene_LB <-         ##### Gene location lower bound
  
### The following information is needed on the outcome
Outcome_Name <-       ##### Outcome Name
Exposure_Name <-       ##### Exposure Name
  
library(TwoSampleMR)
library(dplyr)
library(plyr)
library(stringr)
library(ieugwasr)
library(coloc)
### A access token is required
### Further information about the token system can be found here https://api.opengwas.io/api/ 

### The following outputs will be generated
### A table listing all posterior probability values, no of variants


Gene_ColocSum <- data.frame(matrix(ncol = 12, nrow = 0))
y <-  c("nsnps", "PP.H0", "PP.H1" , "PP.H2" , "PP.H3" , "PP.H4" , "Gene" , "A" , "B",  "CHROM" , "LowB" , "UpB" )
colnames(Gene_ColocSum) <- y

####### Analysis

Exposure_GWAS <-  subset(Exposure_GWAS, Exposure_GWAS$EAF < 1 )
Exposure_GWAS <-  subset(Exposure_GWAS, Exposure_GWAS$P.value > 0 )

Exp_GWAS <- format_data(Exposure_GWAS, 
                             type = "exposure",
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col = "P.value",
                             chr_col = "CHROM",
                             pos_col = "POS",
                             effect_allele_col = "Effect_allele",
                             other_allele_col = "Other_allele",
                             eaf_col = "EAF",
                             samplesize_col = "N")   

IV_List <- subset(Exp_GWAS, Exp_GWAS$pval.exposure < 1e-07 )

#### Clump data using TwoSampleMR package- reselect MR instruments
IV_SNPs <- clump_data(IV_List, clump_kb = 10000,
                      clump_r2 = 0.001,
                      clump_p1 = 5E-08,
                      clump_p2 = 5E-08,
                      pop= "EUR")


#### Performs coloc for 500kb region around each MR IV
##### For quantative output datasets

for (k in 1: nrow(IV_SNPs)) {
  
  Lead_SNP <- IV_SNPs [k,] 
  
  SNP <- paste(Lead_SNP$SNP)
  Chromosome <- as.numeric(paste(Lead_SNP$chr.exposure))
  Position <- as.numeric(paste(Lead_SNP$pos.exposure))
  
  MaxV <- Position + 500000
  LowV <- Position - 500000
  
  Exp_Locus <- subset(Exp_GWAS, CHROM  %in% c(paste(Chromosome)))
  
  Exp_Locus <-  subset(Exp_Locus, POS %in% c(paste(LowV):paste(MaxV) ))
  
Outcome_GWAS <-  subset(Outcome_GWAS, Outcome_GWAS$EAF < 1 )
Outcome_GWAS <-  subset(Outcome_GWAS, Outcome_GWAS$P.value > 0 )

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
                        eaf_col = "EAF",
                        samplesize_col = "N") 

Out_GWAS <-  subset(Out_GWAS, chr.outcome == Chromosome )

matched <- intersect( Out_GWAS$SNP, Exp_Locus$SNP)
all <-  union(Out_GWAS$SNP, Exp_Locus$SNP)
non.matched <- all[!all %in% matched]
Out_SNPS <- Out_GWAS[ which( ! Out_GWAS$SNP 
                             %in% non.matched) , ]

Exp_SNPS <- Exp_Locus[ which( Exp_Locus$SNP 
                             %in% matched) , ]

Exp_SNPS <- Exp_SNPS[order(Exp_SNPS$SNP),]
Out_SNPS <- Out_SNPS[order(Out_SNPS$SNP),]


Harmonised_dat <- harmonise_data(exposure_dat = Exp_SNPS, outcome_dat = Out_SNPS)


Exp_Input <- Harmonised_dat[c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
                              "effect_allele.exposure", "other_allele.exposure",
                              "eaf.exposure", "samplesize.exposure")]

Out_Input <- Harmonised_dat[c("SNP", "beta.outcome", "se.outcome", "pval.outcome",
                              "effect_allele.outcome", "other_allele.outcome",
                              "eaf.outcome", "samplesize.outcome")]

Exp_Input <- Exp_Input %>% group_by(SNP) %>% filter(n() == 1)
Out_Input <- Out_Input %>% group_by(SNP) %>% filter(n() == 1)

Exp_Input <-rename(Exp_Input, c("SNP" = "snp"))
Exp_Input <-rename(Exp_Input, c("beta.exposure" = "beta"))
Exp_Input <-rename(Exp_Input, c("se.exposure" = "SE"))
Exp_Input <-rename(Exp_Input, c("pval.exposure" = "pvalues"))
Exp_Input <-rename(Exp_Input, c("eaf.exposure" = "MAF"))
Exp_Input <-rename(Exp_Input, c("samplesize.exposure" = "N"))

Out_Input <-rename(Out_Input, c("SNP" = "snp"))
Out_Input <-rename(Out_Input, c("beta.outcome" = "beta"))
Out_Input <-rename(Out_Input, c("se.outcome" = "SE"))
Out_Input <-rename(Out_Input, c("pval.outcome" = "pvalues"))
Out_Input <-rename(Out_Input, c("eaf.outcome" = "MAF"))
Out_Input <-rename(Out_Input, c("samplesize.outcome" = "N"))

Exp_Input$varbeta <- Exp_Input$SE^2
Out_Input$varbeta <- Out_Input$SE^2

Exp_Input <-na.omit(Exp_Input)
Out_Input <-na.omit(Out_Input)

matched <- intersect( Exp_Input$snp, Out_Input$snp)
Exp_Input <- Exp_Input[ which( Exp_Input$snp %in% matched) , ]
Out_Input <- Out_Input[ which( Out_Input$snp %in% matched) , ]

Exp_List <- as.list(Exp_Input)
Out_List <- as.list(Out_Input)

Exp_List$type = "quant"
Out_List$type = "quant" ###### If you are using a case control outcome this will need adjusting

cdf <- coloc.abf(Exp_List, Out_List, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)

cdf_results <- cdf$results
cdf_summary <- cdf$summary

names(cdf_summary) <- NULL
cdf_summary <- append(cdf_summary, print( paste(Gene_Name))) ##### Gene
cdf_summary <- append(cdf_summary, print( paste(Exposure_Name))) ### dataset A 
cdf_summary <- append(cdf_summary, print( paste(Outcome_Name))) ####### dataset B
cdf_summary <- append(cdf_summary, print( paste(Chromosome) )) ####### CHROM
cdf_summary <- append(cdf_summary, print( paste(LowV))) ### LowB
cdf_summary <- append(cdf_summary, print( paste(MaxV) )) ####### UpB

cdf_summary<- t(cdf_summary)
cdf_summary <- data.frame(cdf_summary)

colnames(cdf_summary) = c("nsnps", "PP.H0", "PP.H1" , "PP.H2" , "PP.H3" , "PP.H4" , "Gene" , "A" , "B" , "CHROM" , "LowB" , "UpB" )

Gene_ColocSum <- rbind(Gene_ColocSum , cdf_summary)

}

##### For case control output datasets
## in addtion to N (total sample size) an addtional column is required Ncase listing the number of cases in the sample
for (k in 1: nrow(IV_SNPs)) {
  
  Lead_SNP <- IV_SNPs [k,] 
  
  SNP <- paste(Lead_SNP$SNP)
  Chromosome <- as.numeric(paste(Lead_SNP$chr.exposure))
  Position <- as.numeric(paste(Lead_SNP$pos.exposure))
  
  MaxV <- Position + 500000
  LowV <- Position - 500000
  
  Exp_Locus <- subset(Exp_GWAS, CHROM  %in% c(paste(Chromosome)))
  
  Exp_Locus <-  subset(Exp_Locus, POS %in% c(paste(LowV):paste(MaxV) ))
  

Outcome_GWAS <-  subset(Outcome_GWAS, Outcome_GWAS$EAF < 1 )
Outcome_GWAS <-  subset(Outcome_GWAS, Outcome_GWAS$P.value > 0 )

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
                        eaf_col = "EAF",
                        samplesize_col = "N",
                        ncase_col = "Ncase") 

Out_GWAS <-  subset(Out_GWAS, chr.outcome == Chromosome )

matched <- intersect( Out_GWAS$SNP, Exp_Locus$SNP)
all <-  union(Out_GWAS$SNP, Exp_Locus$SNP)
non.matched <- all[!all %in% matched]
Out_SNPS <- Out_GWAS[ which( ! Out_GWAS$SNP 
                             %in% non.matched) , ]

Exp_SNPS <- Exp_GWAS[ which( Exp_Locus$SNP 
                             %in% matched) , ]

Exp_SNPS <- Exp_SNPS[order(Exp_SNPS$SNP),]
Out_SNPS <- Out_SNPS[order(Out_SNPS$SNP),]


Harmonised_dat <- harmonise_data(exposure_dat = Exp_SNPS, outcome_dat = Out_SNPS)


Exp_Input <- Harmonised_dat[c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
                              "effect_allele.exposure", "other_allele.exposure",
                              "eaf.exposure", "samplesize.exposure")]

Out_Input <- Harmonised_dat[c("SNP", "beta.outcome", "se.outcome", "pval.outcome",
                              "effect_allele.outcome", "other_allele.outcome",
                              "eaf.outcome", "samplesize.outcome", "ncase.outcome")]

Exp_Input <- Exp_Input %>% group_by(SNP) %>% filter(n() == 1)
Out_Input <- Out_Input %>% group_by(SNP) %>% filter(n() == 1)

Exp_Input <-rename(Exp_Input, c("SNP" = "snp"))
Exp_Input <-rename(Exp_Input, c("beta.exposure" = "beta"))
Exp_Input <-rename(Exp_Input, c("se.exposure" = "SE"))
Exp_Input <-rename(Exp_Input, c("pval.exposure" = "pvalues"))
Exp_Input <-rename(Exp_Input, c("eaf.exposure" = "MAF"))
Exp_Input <-rename(Exp_Input, c("samplesize.exposure" = "N"))

Out_Input <-rename(Out_Input, c("SNP" = "snp"))
Out_Input <-rename(Out_Input, c("beta.outcome" = "beta"))
Out_Input <-rename(Out_Input, c("se.outcome" = "SE"))
Out_Input <-rename(Out_Input, c("pval.outcome" = "pvalues"))
Out_Input <-rename(Out_Input, c("eaf.outcome" = "MAF"))
Out_Input <-rename(Out_Input, c("samplesize.outcome" = "N"))

Exp_Input$varbeta <- Exp_Input$SE^2
Out_Input$varbeta <- Out_Input$SE^2
Out_Input$s <- Out_Input$ncase.outcome/Out_Input$N

Exp_Input <-na.omit(Exp_Input)
Out_Input <-na.omit(Out_Input)

matched <- intersect( Exp_Input$snp, Out_Input$snp)
Exp_Input <- Exp_Input[ which( Exp_Input$snp %in% matched) , ]
Out_Input <- Out_Input[ which( Out_Input$snp %in% matched) , ]

Exp_List <- as.list(Exp_Input)
Out_List <- as.list(Out_Input)

Exp_List$type = "quant"
Out_List$type = "cc" ###### If you are using a case control outcome 

cdf <- coloc.abf(Exp_List, Out_List, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)

cdf_results <- cdf$results
cdf_summary <- cdf$summary

names(cdf_summary) <- NULL
cdf_summary <- append(cdf_summary, print( paste(Gene_Name))) ##### Gene
cdf_summary <- append(cdf_summary, print( paste(Exposure_Name))) ### dataset A 
cdf_summary <- append(cdf_summary, print( paste(Outcome_Name))) ####### dataset B
cdf_summary <- append(cdf_summary, print( paste(Chromosome) )) ####### CHROM
cdf_summary <- append(cdf_summary, print( paste(LowV))) ### LowB
cdf_summary <- append(cdf_summary, print( paste(MaxV) )) ####### UpB

cdf_summary<- t(cdf_summary)
cdf_summary <- data.frame(cdf_summary)

colnames(cdf_summary) = c("nsnps", "PP.H0", "PP.H1" , "PP.H2" , "PP.H3" , "PP.H4" , "Gene" , "A" , "B" , "CHROM" , "LowB" , "UpB" )

Gene_ColocSum <- rbind(Gene_ColocSum , cdf_summary)

}
