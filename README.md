# VaD_DTMR2025
Taylor-Bateman et al., 2025 - Repurposing drugs for the prevention of vascular dementia: Evidence from drug target Mendelian randomization 

File descriptions- 
Cis_DTMR.R - runs MR for an exposure and outcome using instruments from the chosen cis-acting gene region
CisLT_DTMR.R - Runs cis-MR at a lower IV selection threshold
Cis&Trans_DTMR.R - Runs a MR analysis using both cis and trans IVs
GColoc_DTMR.R - Runs a coloc analysis for the chosen gene region
SNPColoc_DTMR.R - Runs a coloc analysis for +/-500kb around each Cis/Trans MR instrument
Ex_Exposure.txt - A simulated exposure (not actual data) for demonstrative use
Ex_Outcome.txt - A simulated outcome (not actual data) for demonstrative use, for coloc this is a quantitative dataset

## Software used-
R (version 4.4.0)
R packages-
TwoSampleMR (version 0.6.8)
stringr (version 1.5.1)
ieugwasr (version 1.0.3)
coloc (version 5.2.3)
dplyr (version 1.1.4)
plyr (version 1.8.9)

## Instructions for Demo-
Load Exposure and Outcome datasets as Exposure_GWAS and Outcome_GWAS respectively
The following information can be added to the R script in the relevant sections if required

Gene_Name <- "PCSK9"       ##### Gene identifer 

Gene_UB <-  55064852      ##### Gene location upper bound 

Gene_LB <- 55039548         ##### Gene location lower bound

Chromosome <- 1     ##### Gene Chromosome

Outcome_Name <- "VAD_Test"     ##### Outcome Name

Exposure_Name <- "PCSK9"      ##### Exposure Name

Each Script should take less than 10 minutes to run for each exposure/outcome comparison












