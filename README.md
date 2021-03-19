
# AT-HLA 
  


# 1. Introduction

The Human Leukocyte Antigen System (HLA) is a group of proteins encoded by the Major Histocompatibility Complex (MHC) gene complex in humans located in the 6<sup>th</sup> chromosome.

<p align="center">
  <img src="https://www.researchgate.net/profile/Thayne-Sweeten/publication/221914655/figure/fig1/AS:393956899737611@1470938052496/A-simplified-map-of-the-HLA-region-on-human-chromosome-6.png">
</p>
<p align="center">Figure 1: Simplified map of the HLA region on human chromosome 6.</p>

These cell-surface proteins are responsible for the regulation of the immune system and for triggering immune responses when a foreign organism is encountered. There are three different classes of MHC encoded by multiple regions of HLA with different functions:
1. MHC Class I: encoded in HLA-A, B, and C. These proteins are responsible for presenting peptides from inside the cell. For instance, if the cell is infected by a foreign organism (e.g. a virus) the HLA system will bring fragments of the virus to the surface of the cell so that it is identified by the immune system as an infected cell and is destroyed by killer T-cells.
2. MHC Class II: encoded in HLA-DP, DM, DO, DQ and DR. These proteins present antigens from outside the cell to T-cells, and stimulate the multiplication of T-helper cells which in turn stimulate antibody-producing B-cells to produce antibodies to that specific antigen.
3. MHC Class III: encodes components of the complement system. 

The HLA system is highly polymorphic - that is, there are many different types of alleles. Hence, the importance of studying the association of specific genes with a condition.

# 2. Getting Started
## 2.1. Requirements 
AT-HLA runs on:
> R 4.0.2 or higher <br>

with the following packages:
> jsonlite 1.7.0 <br>
> tidyverse 1.3.0 <br>
> readr 1.3.1 <br>
> data.table 1.12.8 <br>
> xlsx 0.6.3 <br>
> zeallot 0.1.0 <br>
> epitools 0.5.10.1 <br>
> plyr 1.8.6 <br>
> haplo.stats 1.8.6 <br>

After installing the aforementioned packages, clone the repository anywhere in your local machine with <br>
```
git clone https://github.com/vicentepese/AT-HLA
``` 

## 2.2. Initilization

Prior to utilizing AT-HLA, the Input Data, and Settings must be initialized. 

### 2.2.1. Input Data

The HLA data (i.e. Input Data) must be provided in a *.csv* file. In addition, it must follow a specific format whereby each row corresponds to a unique instance or subject, and pairs of alleles are separate in contiguous columns. Please consider the following:
- Subjects' unique identifiers must be stored in a column named *<span>sample</span>.id*.
- Each pair of alleles must be named in the format of "*LOCUS*.1", "*LOCUS*.2" whereby *LOCUS* is any of the HLA loci. The total number of columns should be *1+2N* where N is the number of loci &ndash; columns with different names will be ignored.

A typical Input Data should look like:<br>
|<span>sample</span>.id   | A.1     | A.2   | B.1   | B.2   | ..... | DRB1.1 | DRB1.2 |
| ----------------------- | ------- | ----- | ----- | ----- | ----- | -----  | -----  |
| ID001                   | 01:01   | 01:02 | 13:03 | 03:04 | ..... | 07:07  | 08:09  |
| ID002                   | 06:03   | 01:02 | 04:01 | 02:04 | ..... | 01:01  | 02:03  |
| ID003                   | 01:02   | 02:02 | 07:01 | 03:04 | ..... | 09:08  | 09:00  |
| .....                   | .....   | ..... | ..... | ..... | ..... | 13:01  | 13:01  |
| ID001                   | 01:01   | 01:02 | 13:03 | 03:04 | ..... | 09:10  | 04:04  |

**Note**: AT-HLA does not have a digit precision limit.

### 2.2.2. Imputation Probabilities

An imputation probability cut-off is performed by AT-HLA &ndash; every subject below a probability threshold will be excluded. This threshold is computed by locus &ndash; that is, if a subject has an HLA-A imputation probability lower than the threshold, it will be excluded in the analysis of the HLA-A locus. However, if the same subject has an HLA-DRB1 imputation probability higher than the threshold, it will not be excluded from the HLA-DRB1 analysis. <br>
<br>
Identically to the Input Data,  a *.csv* file in which each row corresponds to a unique subject must be provided. Please consider:
- Subjects' unique identifiers must be stored in a column named *<span>sample</span>.id*. Such identifiers must match those of the Input Data.
- Each locus imputation probability must be named "prob.*LOCUS*" whereby *LOCUS* is any of the HLA loci. The total number of columns should be *1+N* where *N* is the number of loci.

A typical Imputation Probabilities file should look like:
|<span>sample</span>.id   | prob.A  | prob.B | ..... | prob.DRB1 |
| ----------------------- | ------- | -----  | ----- | -----  |
| ID001                   | 0.987   | 1.000  | ..... | 0.785  |
| ID002                   | 0.778   | 0.998  | ..... | 0.995  |
| ID003                   | 0.122   | 0.255  | ..... | 0.653  |
| .....                   | .....   | .....  | ..... | .....  |
| ID001                   | 0.236   | 0.356  | ..... | 0.846  |

### 2.2.3. Covariates

A covariates *.csv* file in which each row corresponds to a unique subject must be provided. This file contains the phenotype of the subjects. Please consider:
- Subjects' unique identifiers must be stored in a column named *<span>sample</span>.id*. Such identifiers must match those of the Input Data.
- Phenotype type must be stored in a column named *pheno*, where 1=cases and 0=controls, or 2=cases and 1=controls.
- *Optional*; only for Generalized Linear Models (GLM) analyses: Principal Components must be provided in the form of *PCX* whereby *X* is the Principal Component number. 

A typical Covariates file should look like this
|<span>sample</span>.id   | PC1     | PC2    | ..... | pheno  |
| ----------------------- | ------- | -----  | ----- | -----  |
| ID001                   | 0.005   | -0.069 | ..... | 1      |
| ID002                   | 0.655   | 0.331  | ..... | 1      |
| ID003                   | -0.362  | 0.253  | ..... | 0      |
| .....                   | .....   | .....  | ..... | .....  |
| ID001                   | 0.158   | -0.258 | ..... | 1      |

**NOTE**: At least the first three Principal Components (i.e., PC1, PC2, and PC3) must be included in the file &ndash; otherwise, scripts requiring covariates will not be functional.

### 2.2.4 Ethnicity (optional)

AT-HLA allows ethnic-specific analyses. Ethnicity should be inputted to AT-HLA in the same format as Input Data, whereby each row corresponds to a single subject's instance. Please consider:
- Subjects' unique identifiers must be stored in a column named *<span>sample</span>.id*. Such identifiers must match those of the Input Data.
- Subjects' ethnicity must be stored in a column named *Population*.

A typical ethnicity file should look like:
|<span>sample</span>.id   | Population     | 
| ----------------------- | ------- | 
| ID001                   | EUR     | 
| ID002                   | EUR     | 
| ID003                   | EAS     |
| .....                   | .....   | 
| ID001                   | SAS     | 

Note that AT-HLA does not require specific labels for *Population*. The provided example uses the ethnicities labeled by [1000 Genomes](https://www.internationalgenome.org/category/population/)

### 2.2.3. Settings

AT-HLA follows a *settings* logic. This means that paths, files, and variables are stored in `settings.json`. To use AT-HLA, please fill-up the following sections in `settings.json`:
```
file:
    HLA_Data      [string]        Full path to the HLA data to be analyzed (i.e. Input Data)
    covars        [string]        Full path to the covariates of the HLA data
    probs         [string]        Full path to the imputation probabilities of the HLA data.
    probs_thr     [float]         Probability threshold. See section 2.2.2. If the Input Data is
                                  unimputed, by setting this variable to 0 and no subjects will be excluded.
    freq_thr      [float]         Frequency threshold. Only applied for P-value correction. P-values belonging to alleles for which the carrier/allele frequencies both in cases and controls are lower than the threshold both in cases and controls, will not be corrected. This allows a less stringent P-value correction.
  
```
- *file*:
  - *HLA_Data* (string): Full path to the HLA data to be analyzed (i.e. Input Data)
  - *covars* (string): Full path to the covariates of the HLA data
  - *probs* (string): Full path to the imputation probabilities of the HLA data.
- *prob_thr* (float): Probability threshold. See section 2.2.2. If the Input Data is unimputed, by setting this variable to 0 and no subjects will be excluded.
- *freq_thr* (float): Frequency threshold. Only applied for P-value correction. P-values belonging to alleles for which the carrier/allele frequencies both in cases and controls are lower than the threshold both in cases and controls, will not be corrected. This allows a less stringent P-value correction.

Optional settings:
- *file*: 
  - *ethnicity* (string): Full path to the HLA data subjects' ethnicity.
- *ethnicity* (list of strings): If the ethnicity of subjects is provided, select the specific ethnicities (provided by the ethnicity file) of interest to perform the analysis (e.g.: EUR, SAS, EAS, AFR, AMR).
- *allele2exclude* (list of strings): The alleles written in this list will be removed from the analysis &ndash; that is, only negative subjects will be considered for the analysis. Alleles must be written into the list with a *LOCUS*\*XX:XX nomenclature (e.g. DRB1\*07:01, DRB1\*04:02)

Additional optional settings will be described in each script's section.

# 3. HLA Association Analysis
An HLA association analysis allows identifying specific alleles that may be associated with a determined condition.  

## 3.1 Allele Association Analysis

### 3.1.1 Command:
From directory of the cloned repository: <br>
```
Rscript HLA_Chi2.R
```

### 3.1.2 Description:
Computes a Chi-square and Fisher's Exact Test in alleles and carrier counts. P-values are locus-wise False Discovery Rate (FDR) corrected by the Benjamini–Yekutieli procedure. The counts of homozygous, heterozygous, non-carriers, and frequencies of cases and controls are included in the output. 

### 3.1.3 Outputs and results interpretation:
The script will produce two outputs, a carrier, and an allelic association study:
- Carrier Association Analysis: *Outputs/Chi2/HLA_AnalysisCarriers.xlsx*
- Allele Association Analysis:  *Outputs/Chi2/HLA_AnalysisAlleles.xlsx*


Both files will contain the following common columns:
- *allele*: Allele of study
- *FishersPVAL*: Fisher's Exact Test P-value
- *FishersPVAL_CORR*: Fisher's Exact Test P-value FDR corrected
- *FishersOR*: Fisher's Odd Ratios (OR)
- *FishersLI*: Fisher's OR Lower 95% Confidence Interval
- *FishersUI*: FIsher's OR Upper 95% Confidence Interval
- *ChiPVAL*: Chi-square P-value
- *ChiPVAL_CORR*: Chi-square P-value FDR corrected
- *OR*: Odds Ratio computed directly from the contingency matrix
- *A0 / A1 / A2 in Case/Control*: In order, these columns correspond to the non-carrier, heterozygous, and homozygous counts in cases and controls.
- *carrier/alleleFreq*: Carrier or allele frequency.
- *carrier/alleleCount*: Count of the carrier or allele specified in *allele*.
- *carrier/alleleTotal*: Total number of carriers or alleles included in the analysis.

## 3.2 Logistic Model Analysis

### 3.2.1 Command
From directory of the cloned repository: <br>
```
Rscript HLA_glm.R
```

### 3.2.2 Description 
Fits a Generalized Logistic Model (GLM) to each allele, controlling by the first three Principal Components. P-values are FDR corrected and through BY procedure. Counts of homozygous, heterozygous, non-carriers, and frequencies in cases and controls are computed as well. 

### 3.2.3 Outputs and results interpretation
The script will produce two outputs, a carrier, and an allelic GLM association analysis:
- Carrier GLM Association Study: *Outputs/GLM/HLA_Carriers.xlsx*
- Allele GLM Association Study: *Outputs/GLM/HLA_Alleles.xlsx*

Both files will contain the following common columns:
- *allele*: Allele of study
- *allele.COEF*: GLM coefficient of the allele 
- *Intercept.pval*: P-value of the intercept
- *allele.pval*: P-value of the allele
- *allele.pval.CORR*: P-value of the allele FDR corrected
- *PC1/2/3.pval*: P-alue of the PC1,PC2, and PC3

For allele and carrier counts column reference, please see Section 3.1: Outputs and results interpretation.

### 3.2.4 Additional settings
The following additional settings can be adjusted through `settings.json`:
- *allele2control* (list of strings): Controls in the GLM model for the specified alleles in the list. 

**Note**: The nomenclature *must* be LOCUS\*XX:XX (see Section 2.2.3: *allele2exclude*). P-values of the controlled alleles will be displayed in the output as *LOCUS\*XX:XX.pval*

## 3.3 Iterative Logistic Model Analysis

### 3.3.1 Command
From directory of the cloned repository: <br>
```
Rscript HLA_glm_iter.R
```

### 3.3.2 Description
Computes a Generalized Logistic model for each allele, and iteratively will control for the most BY-FDR corrected significant allele of the previous iteration until no significant alleles are left. Counts of homozygous, heterozygous, non-carriers, and frequencies in cases and controls are computed as well.

### 3.3.3. Outputs and results interpretation
The script will produce a carrier iterative GLM association analysis:
- Iterative Carrier GLM Association: *Outputs/GLM/HLA_GLM_Carriers_iter.xlsx*

The output format is the same as in Section 3.2.3. In addition, the output includes an additional Excel Sheet named *Significant_alleles* that includes the most significant allele for each iteration, from top to bottom.  

### 3.3.4 Additional settings
The following additional settings can be adjusted through `settings.json`:
- *allele2control* (list of strings): Controls in the GLM model for the specified alleles in the list for all iterations (including the first iteration).
- *skipAllele* (list of strings): Ignores the alleles defined in the list. This is because sometimes alleles can be collinear, and previously controlled alleles will become non-significant &ndash; skipping such alleles will avoid redundant results.

**Note**: The nomenclature *must* be LOCUS\*XX:XX (see Section 2.2.3: *allele2exclude*). 

# 4. Haplotype Analysis

## 4.1. Haplotype count 

### 4.1.1. Command 
From directory of the cloned repository: <br>
```
Rscript HLA_haplotype.R
```

### 4.1.2 Decription
Provided an haplotype, it performs a count of heterozygous and homozygous cases and controls, and will test agains a reference of non-carriers.

### 4.1.3 Settings
Please, in addition to the necessary settings described in Section 2.2.3 include:
- *Haplotype* (list of strings): haplotypes, i.e list of alleles.

**Note**: The nomenclature *must* be LOCUS\*XX:XX (see Section 2.2.3: *allele2exclude*). 

### 4.1.4  Outputs and results interpretation
The script will produce a haplotype count association output:
- Haplotype association: *Outputs/Haplotype/HLA_Haplotype.xlsx*

The file contains the following columns:
- *LOCUS.1/.2*: Allele 1 and allele 2 
- *Ncases / Ncontrols*: Number of cases/controls carrying the haplotype
- *FreqCases / FreqControls*: Frequency of cases/controls carrying the haplotype
- *OR*: Odds Ratio computed from directly from the contingency matrix
- *chi.sq*: Chi-square P-value.

## 4.1. Expectation Maximization

### 4.1.1. Command 
From directory of the cloned repository: <br>
```
Rscript HLA_EM.R
```

### 4.1.2 Decription
Provided an set of loci, it performs the Expectation Maximization (EM) algorithm and predict haplotype frequencies.

### 4.1.3 Settings
Please, in addition to the necessary settings described in Section 2.2.3 include:
- *EM_haplo* (list of strings): haplotypes, i.e list of loci.

**Note**: The nomenclature *must* be LOCUS (e.g, DRB1, DPB1)

### 4.1.4  Outputs and results interpretation
The script will produce a EM haplotype count association output:
- EM Haplotype association: *Outputs/Haplotype/HLA_EM.xlsx*

The file contains the following columns:
- *LOCUS*: Allele of the locus/loci 
- *Count.Cases/Controls*: Number of cases/controls carrying the haplotype
- *FreqCases / FreqControls*: Frequency of cases/controls carrying the haplotype
- *OR*: Odds Ratio computed from directly from the contingency matrix
- *Chi2*: Chi-square P-value.
- *RefCases/Controls*: Number of cases and controls not carrying the haplotype.

# 5. Zygosity Analysis

# 6. Amino Acid Association Analysis



