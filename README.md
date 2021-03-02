# AT-HLA: Automatic HLA Analysis


# Table of Contents

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

# 2. Initialization
## 2.1. Requirements 
AT-HLA runs on:
> R 4.0.2 or superior <br>

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

## 2.2 Getting started 

Clone the repository anywhere in your local machine with <br>
```
git clone https://github.com/vicentepese/AT-HLA
``` 
Copy the HLA calls and 
Copy your HLA calls into the `Data` folder. The file *must* be a *.csv* file whereby each row is a subject and the following columns:
- *sample.*<span></span>*id*: Subjects' unique identifier.
- 



