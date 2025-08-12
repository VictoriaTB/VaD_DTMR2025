# VaD_DTMR2025
Taylor-Bateman et al., 2025 - Repurposing drugs for the prevention of vascular dementia: Evidence from drug target Mendelian randomization 

# MR & Colocalization Analysis Scripts

## File Descriptions
- **`Cis_DTMR.R`** – Runs MR for an exposure and outcome using instruments from the chosen cis-acting gene region.  
- **`CisLT_DTMR.R`** – Runs cis-MR at a lower IV selection threshold.  
- **`Cis&Trans_DTMR.R`** – Runs MR using both cis and trans IVs.  
- **`GColoc_DTMR.R`** – Runs a coloc analysis for the chosen gene region.  
- **`SNPColoc_DTMR.R`** – Runs a coloc analysis for ±500kb around each cis/trans MR instrument.  
- **`Ex_Exposure.txt`** – Simulated exposure dataset (demo only, not actual data).  
- **`Ex_Outcome.txt`** – Simulated outcome dataset (demo only, quantitative data for coloc).  
- **`Ex_Outputs.xlsx`** – Simulated results from the demo for comparison.

---

## R Environment and Packages
Analysis performed using:
- **R version**: 4.4.0  
- **RStudio version**: 2025.05.0+496 (“Mariposa Orchid” Release)

### Required R Packages

| Package       | Version |
|---------------|---------|
| TwoSampleMR   | 0.6.8   |
| stringr       | 1.5.1   |
| ieugwasr      | 1.0.3   |
| coloc         | 5.2.3   |
| dplyr         | 1.1.4   |
| plyr          | 1.8.9   |

---

## Installing R and RStudio
1. [Download R from CRAN](https://cran.r-project.org/)  
2. [Download RStudio Desktop from Posit](https://posit.co/download/rstudio-desktop/)  

> Installation typically takes less than 30 minutes in total.

---

## Installing Packages
```r
install.packages(c("TwoSampleMR", "stringr", "ieugwasr", "coloc", "dplyr", "plyr"))
```
- Installs all required packages from CRAN in one step.  
- Typically takes less than 5 minutes per package (may be faster if dependencies are already installed).

---

## Loading Packages
After installation, load the packages at the start of your R script:

```r
library(TwoSampleMR)
library(stringr)
library(ieugwasr)
library(coloc)
library(dplyr)
library(plyr)
```

---

## Access Tokens
To use the **OpenGWAS** database with the `TwoSampleMR` package, you must create a free account and obtain an API access token.

- Detailed instructions are available here:  
  [ieugwasr authentication guide](https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication)  
- API information can be found at:  
  [https://api.opengwas.io/api/](https://api.opengwas.io/api/)

---

## Data Formats

The expected GWAS input format (also documented at the start of each script) is:

| Column name     | Description                                           |
|-----------------|-------------------------------------------------------|
| **ID**          | rsID of the SNP (e.g., `rs123456`).                   |
| **CHROM**       | Chromosome number where the SNP is located.           |
| **POS**         | Base-pair position of the SNP (GRCh38/hg38 build).    |
| **Effect_allele** | Allele for which the effect size is reported.        |
| **Other_allele**  | The alternative allele at the same locus.            |
| **Beta**        | Estimated effect size of the effect allele on trait.  |
| **EAF**         | Effect allele frequency (if available).                |
| **SE**          | Standard error of the effect size estimate.            |
| **P.value**     | P-value for the SNP-trait association.                 |

---

## How to Run the Scripts

1. Open the desired R script in R or RStudio.  
2. Edit the parameters at the top of the script to specify:  
   - Gene name  
   - Gene location (chromosome, upper and lower base-pair bounds)  
   - Outcome name  
   - Exposure name  
3. Run or source the script in the R console.

---

## Output Files

**MR Analysis produces two result tables:**  
- `Output_IV_Info` — information about the instrumental variables used in the analysis.  
- `Output_Results` — MR results including MR-Egger, Weighted Median (if applicable), and the F-statistic.

**Colocalization Analysis produces one result table:**  
- `Gene_ColocSum` — posterior probabilities for colocalization hypotheses and number of SNPs analyzed.

Example input and output files are provided in the repository for reference.

---

## Demo Instructions

Load the example exposure and outcome datasets as follows and add the following parameters to the start of the script:

```r
Exposure_GWAS <- read.table("Ex_Exposure.txt", header = TRUE)
Outcome_GWAS  <- read.table("Ex_Outcome.txt", header = TRUE)

# Define parameters
Gene_Name     <- "PCSK9"       # Gene identifier
Gene_UB       <- 55064852      # Gene upper bound (bp)
Gene_LB       <- 55039548      # Gene lower bound (bp)
Chromosome    <- 1             # Chromosome number
Outcome_Name  <- "VAD_Test"    # Outcome name
Exposure_Name <- "PCSK9"       # Exposure name
```

Each script should run in less than 10 minutes per exposure/outcome comparison.

---










