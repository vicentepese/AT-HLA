## ---------------------------
##
## Name: HLA_haplotype.R
##
## Desciption: Given a dataset of HLA calls and a list of alleles,
##              it computes an haplotype analysis.
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
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(HIBAG)
library(parallel)
library(corrplot)
library(randomForest)
library(xlsx)
library(plyr)
library(epitools)

########### INITIALIZATION ############

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
haplotype <- settings$Haplotype %>% unlist()

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

# Delete files to allow output to be written
file.names <- list.files(settings$Output$Haplotype, full.names = TRUE)
file.names <- file.names[grepl(file.names, pattern = "HLA_Haplotype.")]
file.remove(file.names)

################## COUNT HAPLOTYPES ################## 

countHaplo = function(settings, HLA.df, covars.df){
  
  # Create haplotype variables 
  haplotype <- settings$Haplotype %>% unlist()
  
  # Parse alleles 
  loci <- haplotype %>% lapply(function(x) x %>% strsplit("\\*") %>% unlist() %>% head(n=1)) %>% unlist()
  probs.loci.id <- paste0(rep("prob.", length(loci)), loci)
  alleles.haplo <- haplotype %>% lapply(function(x) x %>% strsplit("\\*") %>% unlist() %>% tail(n=1)) %>% unlist()
  
  # Count number of cases
  covars.df <- covars.df %>% filter(sample.id %in% HLA.df$sample.id)
  NCases <- covars.df$pheno %>% table() %>% .['1']; Ncontrols <- covars.df$pheno %>% table() %>% .['2']
  
  # Initialize for loop
  hetero.haplo <- HLA.df ; homo.haplo <- HLA.df
  
  # For each allele, filter the dataframe
  for (i in 1:length(alleles.haplo)){
    
    # Parse locus and allele 
    locus <- loci[i]; allele <- alleles.haplo[i]
    
    # Filter between carriers and homozygous
    hetero.haplo <- hetero.haplo %>% filter(get(paste0(locus,".1")) == allele | get(paste0(locus, ".2")) == allele)
    homo.haplo <- homo.haplo %>% filter(get(paste0(locus,".1")) == allele & get(paste0(locus, ".2")) == allele)
    
  }
  
  # Remove homozygous from carriers to parse heterosygous
  hetero.haplo <- hetero.haplo %>% filter(sample.id %notin% homo.haplo$sample.id)
  
  # Get reference
  ref.haplo <- HLA.df %>% filter(sample.id %notin% hetero.haplo$sample.id & sample.id %notin% homo.haplo$sample.id)
  
  # Parse alleles to exclude
  allele2control <- settings$allele2control %>% unlist()
  loci <- allele2control %>% lapply(function(x) x %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)) %>% unlist()
  alleles <- allele2control %>% lapply(function(x) x %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)) %>% unlist()
  
  # Filter out reference based on allele 2 exclude in setttings
  for (i in 1:length(loci)){
    locus.id <- c(paste0(loci[i], ".1"), paste0(loci[i], ".2"))
    ref.haplo <- ref.haplo %>% filter(get(locus.id[1]) != alleles[i] & get(locus.id[2]) != alleles[i])
  }
  
  # Parse cases and controls ids
  cases.ids <- covars.df %>% filter(pheno ==1) %>% .["sample.id"] %>% unlist(); Ncases <- cases.ids %>% length()
  controls.ids <- covars.df %>% filter(pheno ==0) %>% .["sample.id"] %>% unlist(); Ncontrols <- controls.ids %>% length()
  
  # Count cases and controls for each of the haplotypes, and the reference 
  NcasesHetero <- hetero.haplo %>% filter(sample.id %in% cases.ids) %>% nrow(); 
  NcontrolsHetero <- hetero.haplo %>% filter(sample.id %in% controls.ids) %>% nrow();
  NcasesHomo <- homo.haplo %>% filter(sample.id %in% cases.ids) %>% nrow();
  NcontrolsHomo<- homo.haplo %>% filter(sample.id %in% controls.ids) %>% nrow();
  NcasesRef <- ref.haplo %>% filter(sample.id %in% cases.ids) %>% nrow();
  NcontrolsRef <- ref.haplo %>% filter(sample.id %in% controls.ids) %>% nrow();
  
  
  # Compute frequencies 
  FreqCasesHetero <- NcasesHetero / NCases *100; FreqControlsHetero <- NcontrolsHetero / Ncontrols *100
  FreqCasesHomo <- NcasesHomo / NCases*100; FreqControlsHomo <- NcontrolsHomo / Ncontrols*100;
  FreqCasesRef <- NcasesRef / NCases*100; FreqControlsRef <- NcontrolsRef / Ncontrols*100;
  
  # Create dataframe
  haplo.df <- data.frame(matrix(ncol = length(alleles.haplo)*2, nrow = 3))
  colnames(haplo.df) = c(paste0(rep(loci, 2), rep(c(".1",".2"), length(loci))))
  haplo.df[1,] <- c(rep(haplotype, 2))
  haplo.df[2,] <- c(haplotype, rep("X", length(haplotype)))
  haplo.df[3,] <- rep("X", ncol(haplo.df))
  
  # Append cases and frequencies
  haplo.df$Ncases <- c(NcasesHomo, NcasesHetero, NcasesRef);
  haplo.df$FreqCases <- c(FreqCasesHomo, FreqCasesHetero, FreqCasesRef)
  haplo.df$Ncontrols <- c(NcontrolsHomo, NcontrolsHetero, NcontrolsRef)
  haplo.df$FreqControls <- c(FreqControlsHomo, FreqControlsHetero, FreqControlsRef)
  
  # Compute Odd-ratios and Chi-square
  chi.sq <- c(); OR <- c()
  for (i in range(1,2)){
    
    # Create contingency table 
     cont.table <- matrix(c(haplo.df[i,"Ncases"], haplo.df[3,"Ncases"], haplo.df[i,"Ncontrols"],  haplo.df[3,"Ncontrols"]), ncol = 2)
     
    # Compute odd ratio and append to vector
     OR <- c(OR,(cont.table[1]*cont.table[4])/(cont.table[2]*cont.table[3]))
     
     # Compute Chi square and append to vector 
     chi.sq <- c(chi.sq, chisq.test(cont.table) %>% .["p.value"])
  }
  
  # Append NAs for reference 
  OR <- c(OR, NA); chi.sq <- c(chi.sq, NA)
  
  # Append to dataframe 
  haplo.df$OR <- OR; haplo.df$chi.sq <- chi.sq
  
  # Return dataframe 
  return(haplo.df)
  
}

################## HAPLOTYPE ANALYSIS ################## 

# Create haplotype variables 
haplotype <- settings$Haplotype %>% unlist()

# Parse alleles 
loci <- haplotype %>% lapply(function(x) x %>% strsplit("\\*") %>% unlist() %>% head(n=1)) %>% unlist()
probs.loci.id <- paste0(rep("prob.", length(loci)), loci)

# Filter out HLA calls wiht low imputation probabilities 
probs.df.filt  <- paste0(rep())
probs.df.filt <- probs.df %>% filter_at(.vars = probs.loci.id, all_vars(.>prob_thr))
HLA.df <- HLA.df %>% filter(sample.id %in% probs.df.filt$sample.id)

# Count haplotypes 
haplo.df <- countHaplo(settings, HLA.df, covars.df)

# Write to excel 
write.xlsx(x = haplo.df, file = paste0(settings$Output$Haplotype, "HLA_Haplotype.xlsx"), col.names = TRUE, row.names = FALSE, append = TRUE)
