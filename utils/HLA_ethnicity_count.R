## ---------------------------
##
## Name: HLA_ethnicity_count.R
##
## Desciption: Given a dataset of HLA calls, an allele and
##    an ethnicity, it counts the number of cases/controls per ethnicity.
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

########### INITIALIZATION ########### 

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
eth.df <- read.csv(settings$file$ethnicity)

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

# Correct pheno for simplicity
if (2 %in% covars.df$pheno %>% table() %>% names()){
  covars.df$pheno <- covars.df$pheno - 1
}

# Append phenotype to HLA calls
HLA.df <- merge(HLA.df, covars.df[,c("sample.id","pheno")], by = "sample.id")

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100
alleles2control <- settings$allele2control

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

########### COUNT PER ETHNICITY ########### 

# Initialize dataframe 
alleleCount.df <- data.frame()

# For each ethnicity, count carriers of given alleles
for(eth in settings$ethnicity %>% unlist()){
  
  # Parse HLA calls and total number of cases and controls
  HLA.df_eth <- HLA.df %>% filter(sample.id %in% eth.df$sample.id[which(eth.df$Population == eth)])
  Ncases <- HLA.df_eth$pheno %>% table() %>% .["1"]; NControls <- HLA.df_eth$pheno %>% table() %>% .["0"]
  
  # Count each of the alleles provided 
  for (allele2control in settings$allele2control %>% unlist()){
    
    # Parse locus and alleles
    L2control <- allele2control %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
    L.ids <- paste0(rep(L2control,2), c(".1",".2"))
    A2control <- allele2control %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)
    
    # Filter low imputation probability 
    probs.df_filt <- probs.df %>% filter(get(paste0("prob.",L2control)) > prob_thr)
    HLA.df_eth <- HLA.df_eth %>% filter(sample.id %in% probs.df_filt$sample.id)
    
    # Get carriers 
    HLA.df_carriers <- HLA.df_eth %>% filter(get(L.ids[1]) == A2control | get(L.ids[2]) == A2control)
    
    # Count cases, controls 
    Acases <- HLA.df_carriers$pheno %>% table() %>% .["1"]; Acontrols <- HLA.df_carriers$pheno %>% table() %>% .["0"]
    
    # Bind to dataframe
    alleleCount.df <- rbind(alleleCount.df, c(Acases, Acontrols, Acases/Ncases*100, Acontrols/NControls*100, Ncases, NControls))
    
    
  }
  
}

# Append ethnicity and write column names 
colnames(alleleCount.df) <- c("Ncases","Ncontrols","FreqCases", "FreqControls", "totalCases","totalControls")
alleleCount.df <- add_column(alleleCount.df, eth = settings$ethnicity %>% unlist(), .before = "Ncases")

# Replace NAs with 0
alleleCount.df[is.na(alleleCount.df)] <- 0

# Write output
write.xlsx(x = alleleCount.df, file = paste0(settings$Output$Utils, "HLA_ethnicity_count.xlsx"), row.names = FALSE)
