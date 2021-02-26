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
library(haplo.stats)

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

################## COUNT HAPLOTYPES ################## 

# Parse loci of interest 
loci <- c("DQA1","DQB1","DRB1", "DRB4")
loci.IDs <- paste0(rep(loci,each = 2), rep(c(".1", ".2"), length(loci)))
HLA.df <- HLA.df %>% merge(covars.df[,c("sample.id", "pheno")])

# Filter out HLA calls wiht low imputation probabilities 
probs.loci.id <- paste0(rep("prob.", length(loci)), loci)
probs.df.filt <- probs.df %>% filter_at(.vars = probs.loci.id, all_vars(.>prob_thr))
HLA.df <- HLA.df %>% filter(sample.id %in% probs.df.filt$sample.id)

# Compute expectation maximization
haplo.em.df <- haplo.em(HLA.df[,loci.IDs], locus.label = loci)

# Get haplotypes and haplotype indexes
haplotypes <- haplo.em.df$haplotype
subj.idx <- haplo.em.df$subj.id
hap1.code <- haplo.em.df$hap1code
hap2.code <- haplo.em.df$hap2code

# Create dataframe with haplotypes 
haplo.df <- cbind(data.frame(sample.id = HLA.df$sample.id[haplo.em.df$subj.id]), haplotypes[hap1.code,])
haplo.df2 <- cbind(data.frame(sample.id = HLA.df$sample.id[haplo.em.df$subj.id]), haplotypes[hap2.code,])
haplo.df <- rbind(haplo.df, haplo.df2)

# Remove duplicates (include only carriers)
haplo.df <- haplo.df[!duplicated(haplo.df),]

# Add pheno 
haplo.df <- merge(haplo.df, covars.df[,c("sample.id", "pheno")], by = "sample.id")

# Count haplotypes 
haplo.count.cases <- haplo.df %>% filter(pheno ==1) %>% count(loci)
colnames(haplo.count.cases)[ncol(haplo.count.cases)] <- "Count.Cases"
haplo.count.controls <- haplo.df %>% filter(pheno == 0) %>% count(loci)
colnames(haplo.count.controls)[ncol(haplo.count.controls)] <- "Count.Controls"
haplo.count <- merge(haplo.count.cases, haplo.count.controls, by = loci, all = TRUE)

# Get reference 
HLA.df.ref <- HLA.df %>% filter(DRB1.1 != "07:01" & DRB1.2 != "07:01")

# Compute Chi-square 



###############3 REFERENCE ###############

tst <- HLA.df %>% filter( (DQA1.1 == "02:01" | DQA1.2 == "02:01") &
                            (DQB1.1 == "02:02" | DQB1.2 == "02:02") &
                            (DRB1.1 == "07:01" | DRB1.2 == "07:01") &
                            (DRB4.1 == "01:01" | DRB4.2 == "01:01") )
