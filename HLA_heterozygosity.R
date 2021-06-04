## ---------------------------
##
## Name: HLA_heterozygosity.R
##
## Desciption: Given a dataset of HLA calls, computes a  
##    heterozygosity study of each allele.
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

########### INITIALIZATION ########### 

# Set working directory
setwd(getSrcDirectory()[1])

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls, covariates 
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Correct pheno for simplicity
if (2 %in% covars.df$pheno %>% table() %>% names()){
  covars.df$pheno <- covars.df$pheno - 1
}

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

# Merge phenotype to HLA calls
HLA.df <- merge(HLA.df, covars.df[,c("sample.id", "pheno")], by = "sample.id")

# Parse HLA calls based on ethnicity, if provided
if (!settings$ethnicity %>% is_empty()){
  
  # Parse IDs in ethnicity/ies
  ethnicity.df <- read.csv(settings$file$ethnicity)
  ethnicity.df.filt <- ethnicity.df %>% filter(Population %in% settings$ethnicity %>% unlist())
  HLA.df <- HLA.df %>% 
    filter(sample.id %in% ethnicity.df.filt$sample.id)
}

# Delete files to allow output to be written
file.names <- list.files(settings$Output$Zygosity, full.names = TRUE)
file.names <- file.names[grepl(file.names, pattern = "HLA_Heterozygosity.xlsx")]
file.remove(file.names)

########### HETEROZYGOSITY ANALYSIS ########### 

# Parse allele to study 
locus <- settings$allele2control %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
locus.ids <- paste0(rep(locus,2), c(".1",".2"))
allele <- settings$allele2control %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)

# Filter HLA calls for low imputation probability 
probs.df_filt <- probs.df %>% filter(get(paste0("prob.",locus)) > prob_thr)
HLA.df <- HLA.df %>% filter(sample.id %in% probs.df_filt$sample.id)

# Parse data with heterozygous allele, and the other alleles
HLA.df <- HLA.df[xor(x = HLA.df[locus.ids[1]] == allele, y = HLA.df[locus.ids[2]] == allele),]

# Get other alleles 
allele.ids <- HLA.df[,c(locus.ids)] %>% as.matrix() %>% as.vector() %>% unique()
allele.ids <- allele.ids[-which(allele.ids == allele)]

# Separate in cases and controls 
HLA.df.cases <- HLA.df %>% filter(pheno == 1)
HLA.df.controls <- HLA.df %>% filter(pheno == 0)

# For eacch allele, compute Chi2
HeteroA.count <- data.frame()
for (A in allele.ids){
  
 # Count carriers in cases and controls, 
  Ncases <- HLA.df.cases %>% filter(get(locus.ids[1]) == A | get(locus.ids[2]) == A) %>% nrow()
  Ncontrols <- HLA.df.controls %>% filter(get(locus.ids[1]) == A | get(locus.ids[2]) == A) %>% nrow()
  
  # Compute contingency table 
  cont.table <- matrix(c(Ncases, nrow(HLA.df.cases)-Ncases, Ncontrols, nrow(HLA.df.controls)-Ncontrols) , nrow = 2)
  
  # Comput chi2 and odds-ratio
  chi2.res <- chisq.test(cont.table) %>% .['p.value']
  OR.res <- oddsratio.wald(cont.table) %>% .[["measure"]]
  OR <- OR.res[2,1]; OR.LI <- OR.res[2,2]; OR.UL <- OR.res[2,3]
    
  # Apped to data.frame 
  HeteroA.count <- rbind(HeteroA.count,
                         c(Ncases, Ncontrols, Ncases/nrow(HLA.df.cases)*100,
                           Ncontrols/nrow(HLA.df.controls)*100, chi2.res,
                           OR, OR.LI, OR.UL) %>% unlist())
}

# Add column names 
colnames(HeteroA.count) <- c("NCases", "Ncontrols","Freq.Cases", "Freq.Controls", "Chi2.pval", "OR", "OR.UpperLim", "OR.LowerLim")

# Append alleles, and adjust p-value
HeteroA.count <- add_column(HeteroA.count, allele=paste0(rep(locus,nrow(HeteroA.count)), rep("*", nrow(HeteroA.count)), allele.ids), .before = "NCases")
HeteroA.count <- add_column(HeteroA.count, Chi2.pval.CORR = p.adjust(HeteroA.count$Chi2.pval, method = "BY"), .after = "Chi2.pval")

# Write output
write.xlsx(HeteroA.count, file = paste0(settings$Output$Zygosity, "HLA_Heterozygosity.xlsx"), row.names = FALSE)

