## ---------------------------
##
## Name: HLA_HardyWeinberg.R
##
## Desciption: Given a dataset of HLA calls, computes a  
##    Hardy Weinberg equilibrium statistical testing.
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
library(epitools)
library(tidyverse)
library(readr)
library(data.table)
library(xlsx)
library(plyr)

########### INITIALIZATION ########### 

# Import settings
settings <- jsonlite::read_json("settings.json")
options(stringsAsFactors = F)

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
file.names <- file.names[grepl(file.names, pattern = "HLA_HardyWeinberg.xlsx")]
file.remove(file.names)


######### ALLELE / CARRIER FREQUENCY #######

computeACFREQ = function(data, locus, Dx){
  
  # Allele frequency 
  A1 <- locus %>% paste('.1', sep = ''); A2 <- locus %>% paste('.2', sep = '')
  alleles <- list(data[, A1], data[, A2]) %>% unlist() 
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  ACFREQ.df <- data.frame(allele = alleles.count %>% names(), 
                          alleleCount = alleles.count %>% as.vector(),
                          alleleFreq = alleles.freq %>% as.vector(), 
                          alleleTotal = nrow(data))
  
  # Remove allele of interest
  ACFREQ.df <- ACFREQ.df %>% filter(allele != settings$allele2control %>% unlist() %>% strsplit('\\*') %>% unlist() %>% .[2])
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',3), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',3), sep = ''))
            return(ACFREQ.df)
          }
  )
}

######## CHI2 / FISHER EXACT FUNCTION #########

computeChi2 = function(ACFREQ.df){
  
  # Initialize 
  chi2.pval.alleles <- c() ; chi2.pval.carriers <- c();
  chi2.OR.alleles <- c(); chi2.OR.carriers <- c();
  fishers.pval.alleles <- c(); fishers.pval.carriers <- c()
  fishers.OR.alleles <- c(); fishers.OR.carriers <- c()
  fishers.UI.alleles <- c(); fishers.UI.carriers <- c()
  fishers.LI.alleles <- c(); fishers.LI.carriers <- c(); OR <- c();
  OR.low <-c(); OR.high <-c();
  
  # For each allele in the locus 
  for (A in ACFREQ.df$allele){
    
    ## Allele frequency 
    # Create contingency table  and run chi2 test
    allele.data <- ACFREQ.df %>% filter(allele == A)
    cont.table.allele <- matrix(c(allele.data$alleleCountCase,
                                  allele.data$alleleCountControl,
                                  allele.data$alleleTotalCase - allele.data$alleleCountCase ,
                                  allele.data$alleleTotalControl - allele.data$alleleCountControl - alleleCount),
                                nrow = 2)
    # Tests
    chi2.allele.res <- chisq.test(x = cont.table.allele)
    fishers.allele.res <- fisher.test(cont.table.allele)
    # Chi2
    chi2.pval.alleles <- c(chi2.pval.alleles, chi2.allele.res$p.value) ;
    chi2.OR.alleles <- c(chi2.OR.alleles, chi2.allele.res$statistic %>% unname())
    # Fishers' exact test 
    fishers.pval.alleles <- c(fishers.pval.alleles, fishers.allele.res$p.value);
    fishers.OR.alleles <- c(fishers.OR.alleles, fishers.allele.res$estimate %>% unname())
    fishers.LI.alleles <- c(fishers.LI.alleles, fishers.allele.res$conf.int[1])
    fishers.UI.alleles <- c(fishers.UI.alleles, fishers.allele.res$conf.int[2])
    #OR
    OR.res <- oddsratio.wald(cont.table.allele)
    OR <- c(OR, OR.res$measure[2,1]); 
    OR.low <- c(OR.low, OR.res$measure[2,2]); 
    OR.high <- c(OR.high,OR.res$measure[2,3]);
    
    
  }
  stats.data <- data.frame(allele = ACFREQ.df$allele, 
                           FishersAllelePVAL = fishers.pval.alleles, FishersAlleleOR = fishers.OR.alleles, 
                           FishersAlleleLI = fishers.LI.alleles, FishersAlleleUI = fishers.UI.alleles,
                           ChiAllelePVAL = chi2.pval.alleles, ChiAlleleVAL = chi2.OR.alleles, OR = OR,
                           OR.low = OR.low, OR.high = OR.high)
  stats.data[is.na(stats.data)] <- 0
  
  return(stats.data)
  
}

########### HETEROZYGOSITY ANALYSIS ########### 

# Parse allele to study 
locus <- settings$allele2control %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
locus.ids <- paste0(rep(locus,2), c(".1",".2"))
allele <- settings$allele2control %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)

# Filter HLA calls for low imputation probability 
probs.df_filt <- probs.df %>% filter(get(paste0("prob.",locus)) > prob_thr)
HLA.df <- HLA.df %>% filter(sample.id %in% probs.df_filt$sample.id)

# Parse allele in locus
alleles <- list(HLA.df[,locus.ids[1]], HLA.df[,locus.ids[2]]) %>% unlist() %>% unique()
alleles <- alleles[-which(alleles == allele)]

# Parse heterozyguous cases
HLA.df.cases <- HLA.df %>% filter(pheno == 1); 
HLA.df.cases <- HLA.df.cases[xor(HLA.df.cases[,locus.ids[1]] == allele, HLA.df.cases[,locus.ids[2]] == allele),c(locus.ids)]

# Parse controls and count the number of alleles to be removed in controls
HLA.df.controls <- HLA.df %>% filter(pheno == 0); totalControls <- HLA.df.controls %>% nrow()
alleleCount <-0
for (A in settings$allele2exclude){
  locus <- A %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
  locus.ids <- paste0(rep(locus,2), c(".1",".2"))
  allele <- A %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)
  alleleCount <- alleleCount + table(HLA.df.controls[,c(locus.ids[1])])[allele] %>% unname() +
    table(HLA.df.controls[,c(locus.ids[2])])[allele] %>% unname() 
}


# Count alleles in cases and controls 
ACFREQ.cases <- computeACFREQ(HLA.df.cases, locus, "case"); totalCases <-unique(ACFREQ.cases$alleleTotalCase)
ACFREQ.controls <- computeACFREQ(HLA.df.controls, locus, "control");  totalControls <- unique(ACFREQ.controls$alleleTotalControl)*2

# Merge 
ACFREQ.df <- merge(ACFREQ.cases, ACFREQ.controls, by = "allele", all = TRUE); ACFREQ.df[is.na(ACFREQ.df)] <- 0
ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df)); ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df))

# Exclude alleles from analysis to avoid stringent P-value correction
for (A in settings$allele2exclude){
  locus <- A %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
  locus.ids <- paste0(rep(locus,2), c(".1",".2"))
  allele.exc <- A %>% unlist() %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)
  ACFREQ.df <- ACFREQ.df %>% filter(allele != allele.exc)
}

# Compute statistical test
stats.data <- computeChi2(ACFREQ.df = ACFREQ.df)

# Correct p-values
stats.data <- add_column(stats.data, FisherAllelePVAL.CORR = p.adjust(stats.data$FishersAllelePVAL), .after = "FishersAllelePVAL")
stats.data <- add_column(stats.data, ChiAllelePVAL.CORR = p.adjust(stats.data$ChiAllelePVAL), .after = "ChiAllelePVAL")

# Write output
HLA.alleles.df <- merge(stats.data[,c(1,which(grepl('Allele', colnames(stats.data))),which(grepl("OR",colnames(stats.data))))],
                        ACFREQ.df[,c(which(grepl(paste(c('allele'), collapse = '|'), colnames(ACFREQ.df))))],
                        by = 'allele')
# Write output
write.xlsx(x = HLA.alleles.df, file = paste0(settings$Output$Zygosity, "HLA_HardyWeinberg.xlsx"), row.names = FALSE)
