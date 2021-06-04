## ---------------------------
##
## Name: HLA_zygosity.R
##
## Desciption: Given a dataset of HLA calls, computes a  
##    zygosity study of each allele.
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
file.names <- file.names[grepl(file.names, pattern = "HLA_Zygosity.xlsx")]
file.remove(file.names)

########### COMPUTE ALLELE/CARRIER COUNT/FREQUENCIES ########### 

computeACFREQ = function(data, locus, Dx){

  # Parse alleles
  A1 <- locus %>% paste0('.1'); A2 <- locus %>% paste0('.2')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()

  # Count carrier counts and freqencies
  carriers <- data[,c(A1, A2)]
  carriers.levels <- list(data[, A1], 
                          data[, A2]) %>% unlist() %>% levels()
  carriers.unique <- apply(carriers, 1, function(x) unique(x)) %>% unlist() %>% as.factor()
  carriers.count <- table (carriers.unique); carriers.count[c(carriers.levels %>% setdiff(carriers.count %>% names()))] <- 0
  carriers.freq <- carriers.count /nrow(data) * 100
  carrier.df <- data.frame(allele = carriers.count %>% names(), 
                           carrierCount = carriers.count %>% as.vector(),
                           carrierFreq = carriers.freq %>% as.vector(),
                           carrierTotal = nrow(data))
  
  # Count heterozygous, homozygous and non-carrier for each allele
  A0 <- c(); A1 <- c(); A2 <- c();
  for (A in levels(as.factor(alleles))){
    
    # Parse data and count
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    
    # Count Hetero, Homo. and non-missing
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  
  # Replace NAs by 0s
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(as.factor(alleles)), A0 = A0, A1 = A1, A2 = A2)
  
  # Compute frequencies 
  HH.data$FreqHetero <- HH.data$A1 / (HH.data$A0 + HH.data$A1 + HH.data$A2) *100
  HH.data$FreqHomo <- HH.data$A2/ (HH.data$A0 + HH.data$A1 + HH.data$A2)*100
  
  # Change column names and return based on diagnosis
  switch (Dx,
          'case' = {
            colnames(HH.data) <- c('allele', paste0(colnames(HH.data) %>% tail(n=-1), rep("Case", length(colnames(HH.data))-1)))
            return(HH.data)
          },
          'control' = {
            colnames(HH.data) <-c('allele', paste0(colnames(HH.data) %>% tail(n=-1), rep("Control", length(colnames(HH.data))-1)))
            return(HH.data)
          }
  )
}


########### COMPUTE CHI SQUARE ########### 

computeChi2 = function(ACFREQ.df){
  
  # Initialize for loop
  ACFREQ.df$Chi2.Hetero.pval <- NA; ACFREQ.df$OR.Hetero <- NA; ACFREQ.df$Chi2.Homo.pval <- NA; ACFREQ.df$OR.Homo <- NA
  
  # Compute Chi2 and OR for each allele
  for (i in 1:nrow(ACFREQ.df)){
    
    # Heterozygous
    cont.table.hetero <- matrix(c(ACFREQ.df$A1Case[i], ACFREQ.df$A1Control[i],
                                  ACFREQ.df$A0Case[i], ACFREQ.df$A0Control[i]), nrow = 2)
    ACFREQ.df$Chi2.Hetero.pval[i] <-  chisq.test(cont.table.hetero) %>%.["p.value"] %>% unlist() %>% unname()
    ACFREQ.df$OR.Hetero[i] <- (cont.table.hetero[1]*cont.table.hetero[4])/(cont.table.hetero[2]*cont.table.hetero[3])
    
    # Homozygous
    cont.table.homo <- matrix(c(ACFREQ.df$A2Case[i], ACFREQ.df$A2Control[i],
                                ACFREQ.df$A0Case[i], ACFREQ.df$A0Control[i]), nrow = 2)
    ACFREQ.df$Chi2.Homo.pval[i] <- chisq.test(cont.table.homo) %>%.["p.value"] %>% unlist() %>% unname()
    ACFREQ.df$OR.Homo[i] <- (cont.table.homo[1]*cont.table.homo[4])/(cont.table.homo[2]*cont.table.homo[3])
  }
  
  # Return dataframe with ORs and Chi2s p-values
  return(ACFREQ.df)
  
}


########### ZYGOSITY ANALYSIS ########### 

# Get cases and controls and separate datasets
cases.ids <- covars.df$sample.id[which(covars.df$pheno ==1)]
controls.ids <- covars.df$sample.id[which(covars.df$pheno ==0)]
data.cases <- HLA.df %>% filter(sample.id %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id %in% controls.ids)

# For each locus 
loci <- colnames(HLA.df)[grepl(colnames(HLA.df), pattern = "\\.1")]
loci <- loci %>% lapply(function(x) strsplit(x, split = "\\.") %>% unlist() %>% head(n=1)) %>% unlist() 
idx <- 1
for (locus in loci){
  
  # Verbose
  print(paste0("Current locus: ", locus))
  
  # Filter out subjects with imputation probability threshold
  probs.df_filt <- probs.df %>% filter(get(paste0("prob.", locus)) > prob_thr)
  data.cases.filt <- data.cases %>% filter(sample.id %in% probs.df_filt$sample.id)
  data.controls.filt <- data.controls %>% filter(sample.id %in% probs.df_filt$sample.id)
  
  # Compute allele and carrier counts and frequencies
  ACFREQ.cases <- computeACFREQ(data.cases.filt, locus, 'case');
  ACFREQ.controls <- computeACFREQ(data.controls.filt, locus, 'control');

  # Merge and clean 
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df <- ACFREQ.df %>% filter(allele %notin% c('','00:00'))
  
  # Compute Chi2 and OR
  ACFREQ.df <- computeChi2(ACFREQ.df)
  ACFREQ.df[ACFREQ.df==Inf] = 'Inf'
  
  # Write 
  write.xlsx(x = ACFREQ.df, file = paste0(settings$Output$Zygosity, 'HLA_Zygosity','.xlsx', sep = ''), sheetName = locus,
             col.names = TRUE, row.names = FALSE, append = TRUE)
}

