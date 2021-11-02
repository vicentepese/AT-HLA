## ---------------------------
##
## Name: HLA_zygosity.R
##
## Desciption: Given a dataset of HLA calls it performs 
##              homo/heterozygosity study of a provided allele.
##
## Author: Vicente Peris Sempere
##
## Year: 2021
##
## Copyright (c) Vicente Peris Sempere, 2021
## Email: vipese@stanford.edu
##
## ---------------------------

# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(xlsx)
library(plyr)
library(epitools)
library(haplo.stats)

############## INITIALIZATION ############## 

# Set working directory
setwd("~/Documents/HLA_association_pipeline")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls, covariates 
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Correct pheno for logistic regression 
if (2 %in% covars.df$pheno %>% table() %>% names()){
  covars.df$pheno <- covars.df$pheno - 1
}

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

# Parse HLA calls based on ethnicity, if provided
if (!settings$ethnicity %>% is_empty()){
  
  # Parse IDs in ethnicity/ies
  ethnicity.df <- read.csv(settings$file$ethnicity)
  ethnicity.df.filt <- ethnicity.df %>% filter(Population %in% settings$ethnicity %>% unlist())
  HLA.df <- HLA.df %>% 
    filter(sample.id %in% ethnicity.df.filt$sample.id)
}

# Parse allele of interest / allele to control
allele2control <- settings$allele2control %>% unlist()

# Ensuer only one allele is present 
if(length(allele2control)>1){
  stop("Only one allele must be provided")
}


############## ZYGOSITY ANALYSIS ############## 

# Append phenotype to HLA calls
HLA.df <- merge(HLA.df, covars.df[,c("sample.id","pheno")], by = "sample.id")

# Parse data with allele 
locus <- allele2control %>% strsplit("\\*") %>% unlist() %>% head(n=1)
locus.id <- paste0(rep(locus,2), c(".1",".2"))
allele <- allele2control %>% strsplit("\\*") %>% unlist() %>% tail(n=1)

# Filter out low imputation probabilities
probs.locus <- paste0("prob.", locus)
probs.df.filt <- probs.df %>% filter(get(probs.locus) > prob_thr)
HLA.df <- HLA.df %>% filter(sample.id %in% probs.df.filt$sample.id)

# Parse heterozygous and homozygous
HLA.subset <- HLA.df %>% filter(get(locus.id[1]) == allele | get(locus.id[2]) == allele)
HLA.homo <- HLA.df %>% filter(get(locus.id[1]) == allele & get(locus.id[2]) == allele)
HLA.hetero <- HLA.subset %>% filter(sample.id %notin% HLA.homo$sample.id)

# Parse population which are allele negative 
HLA.ref <- HLA.df %>% filter(sample.id  %notin% HLA.subset$sample.id)

# Compute Counts and frequencies
zyg.df <- data.frame(c(allele,allele, "X"), c(allele, "X", "X")); colnames(zyg.df) <- locus.id
zyg.df$Ncases <- c(HLA.homo$pheno %>% table() %>% .["1"], HLA.hetero$pheno %>% table() %>% .["1"], HLA.ref$pheno %>% table() %>% .["1"])
zyg.df$FreqCases <- c(zyg.df$Ncases/ HLA.df$pheno %>% table() %>% .["1"]) *100
zyg.df$Ncontrols <- c(HLA.homo$pheno %>% table() %>% .["0"], HLA.hetero$pheno %>% table() %>% .["0"], HLA.ref$pheno %>% table() %>% .["0"])
zyg.df$FreqControls <- c(zyg.df$Ncontrols/ HLA.df$pheno %>% table() %>% .["0"]) *100

# Compute OR
zyg.df$OR <- NA; zyg.df$Chi2 <- NA
for (i in range(1,2)){
  cont.table <- matrix(c(zyg.df$Ncases[i], zyg.df$Ncases[3], zyg.df$Ncontrols[i], zyg.df$Ncontrols[3]), nrow=2)
  chi2.res <- chisq.test(cont.table)
  OR.res <- (cont.table[1]*cont.table[4])/(cont.table[2]*cont.table[3])
  zyg.df$OR[i] <- OR.res
  zyg.df$Chi2[i] <- chi2.res$p.value
}

# Write output
write.xlsx(x = zyg.df, file = paste0(settings$Output$Zygosity, "HLA_zygosity_", allele2control,".xlsx"), row.names = F)
  
