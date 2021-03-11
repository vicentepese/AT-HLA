
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

Prior to utilizing AT-HLA, the Input Data, and Settings must be intialized. 

### 2.2.1. Input Data

The HLA data (i.e. Input Data) must be provided in a *.csv* file. In addition, it must follow a specific format whereby each row correspond to a unique instance or subject, and pairs of alleles are separate in contiguous columns. Please consider the following:
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

An imputation probability cut-off is performed by AT-HLA &ndash; every subject below a probability threshold will be excluded. This threshold is computed by locus - that is, if a subject has an HLA-A imputation probability lower than the threshold, it will be excluded in the analysis of the HLA-A locus. However, if the same subject has a HLA-DRB1 imputation probability higher than the threshold, it will not be excluded of the HLA-DRB1 analysis. <br>
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
- *Optional*; only for Generalized Linear Models (GLM) anylises: Principal Components must be provided in the form of *PCX* whereby *X* is the Principal Component number. 

A typical Covariates file should look like this
|<span>sample</span>.id   | PC1     | PC2    | ..... | pheno  |
| ----------------------- | ------- | -----  | ----- | -----  |
| ID001                   | 0.005   | -0.069 | ..... | 1      |
| ID002                   | 0.655   | 0.331  | ..... | 1      |
| ID003                   | -0.362  | 0.253  | ..... | 0      |
| .....                   | .....   | .....  | ..... | .....  |
| ID001                   | 0.158   | -0.258 | ..... | 1      |

### 2.2.4 Ethnicity (optional)

AT-HLA allows ethnic-specific analyses. Ethnicity should be inputed to AT-HLA in the same format as Input Data, whereby each row corresponds a single subject's instance. Please consider:
- Subjects' unique identifiers must be stored in a column named *<span>sample</span>.id*. Such identifiers must match those of the Input Data.
- Subjects' ethnicity must be stored in a column named *Population*.

A typical ethnicity file should look like:
A typical Covariates file should look like this
|<span>sample</span>.id   | Population     | 
| ----------------------- | ------- | 
| ID001                   | EUR     | 
| ID002                   | EUR     | 
| ID003                   | EAS     |
| .....                   | .....   | 
| ID001                   | SAS     | 

Note that AT-HLA does not require specific labels for *Population*.

### 2.2.3. Settings

AT-HLA follows a *settings* logic. This means that paths, files and variables are stored into `settings.json`. In order to use AT-HLA, please fill-up the following sections in `settings.json`:
- *file*:
  - *HLA_Data* (string): Full path to the HLA data to be analyzed (i.e. Input Data)
  - *covars* (string): Full path to the covariates of teh HLA data
  - *probs* (string): Full path to the imputation probabilities of the HLA data.
- *prob_thr*: Probability threshold. See section 2.2.2. If the Input Data is unimputed, by set this variableto 0 and no subjects will be excluded.
- *freq_thr*: Frequency threshold. Only applied for p-value correction. P-values belonging to alleles for which the carrier/allele frequenci is lower than the threshold both in cases and controls, will not be corrected  (include FDR) *****. This allows a less stringent P-value correction.

Optional settings:
- *file*: 
  - *ethnicity* (string): Full path to the HLA data subjects' ethnicity.
- *ethnicity*: If the ethnicity of subject is provided, select the specific ethnicities (provided by the ethnicity file) of interest to perform the analysis (e.g.: EUR, SAS, EAS, AFR, AMR).

Additional optional settings will be described in each script's section.



