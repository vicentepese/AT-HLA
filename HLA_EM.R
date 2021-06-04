## ---------------------------
##
## Name: HLA_EM.R
##
## Desciption: Given a dataset of HLA calls it performs 
##              expectation maximization to find haplotypes.
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

########### INITIALIZATION ########### 

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

# Filter out the alleles to exclude 
alleles2control <- settings$allele2exclude %>% unlist()
for (allele in alleles2control){
  
  # Parse locus and allele
  locus <- allele %>% strsplit("\\*") %>% unlist() %>% head(n=1)
  A <- allele %>% strsplit("\\*") %>% unlist() %>% tail(n=1)
  
  # Filter HLA calls 
  HLA.df <- HLA.df[which(HLA.df[,paste0(locus,".1")] != A & HLA.df[,paste0(locus,".2")] != A),]
  
}

################# MERGE DRBs ################# 

merge_DRB = function(settings, HLA.df){
  
  # Check that DRB3 4 or 5 are present in the dataset
  if (grepl(colnames(HLA.df), pattern = paste(c("DRB3", "DRB4", "DRB5"), collapse = "|")) %>% any()){
    
    # Get columns with the alleles
    DRB.cols <- HLA.df[, grepl(colnames(HLA.df), pattern = paste(c("DRB3", "DRB4", "DRB5"), collapse = "|"))]
    
    # Create columns on the HLA calls
    HLA.df$DRB.1 <- NA; HLA.df$DRB.2 <- NA
    
    # For each triplet/duplet of alleles
    for (i in 1:nrow(DRB.cols)){
      
      # Get non-missing alleles
      alleles.nonmiss <- which(DRB.cols[i,] != "00:00")
      
      # Get alleles
      locus.ids <- lapply(colnames(DRB.cols)[alleles.nonmiss], function(x) x %>% strsplit("\\.") %>% unlist() %>% head(n=1))%>% unlist()
      As <- paste(locus.ids, DRB.cols[i,alleles.nonmiss], sep = "*")
      
      # If none, or one allele present, correct for format
      if(length(As) == 0){
        As <- c("DRB*00:00", "DRB*00:00")
      } else if (length(As) == 1){
        As <- c(As, "DRB*00:00")
      } 
      
      # Append 
      HLA.df$DRB.1[i] <- As[1]; HLA.df$DRB.2[i] <- As[2]
      
    }
  }
  
  # Return dataframe with DRB3-4-5 merged to DRB
  return(HLA.df)
  
}

###########  EXPECTATION MAXIMIZATION ########### 

haplo_EM = function(settings, HLA.df, probs.df){
  
  # Parse loci of interest and append phenotypes to HLA calls
  loci <- settings$EM_haplo %>% unlist()
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
  subj.ids.idx <- HLA.df$sample.id[subj.idx]
  hap1.code <- haplo.em.df$hap1code
  hap2.code <- haplo.em.df$hap2code
  
  # Create dataframe with haplotypes 
  haplo.df <- cbind(data.frame(sample.id = HLA.df$sample.id[subj.idx]), haplotypes[hap1.code,])
  haplo.df2 <- cbind(data.frame(sample.id = HLA.df$sample.id[subj.idx]), haplotypes[hap2.code,])
  haplo.df <- rbind(haplo.df, haplo.df2)

  # Remove duplicates (include only carriers / extended haplotype)
  haplo.df <- haplo.df[!duplicated(haplo.df),]
  
  # Add phenotype
  haplo.df <- merge(haplo.df, covars.df[,c("sample.id", "pheno")], by = "sample.id")
  
  # Count haplotypes from cases and controls
  haplo.count.cases <- haplo.df %>% filter(pheno ==1) %>% count(loci)
  colnames(haplo.count.cases)[ncol(haplo.count.cases)] <- "Count.Cases"
  haplo.count.controls <- haplo.df %>% filter(pheno == 0) %>% count(loci)
  colnames(haplo.count.controls)[ncol(haplo.count.controls)] <- "Count.Controls"
  haplo.count <- merge(haplo.count.cases, haplo.count.controls, by = loci, all = TRUE)
  
  # Replace NAs with 0s
  haplo.count[is.na(haplo.count)] <- 0
  
  # Compute frequencies 
  Ncases <- HLA.df$pheno %>% table() %>% .["1"]; Ncontrols <- HLA.df$pheno %>% table() %>% .["0"]
  haplo.count <- haplo.count %>% add_column(Freq.Cases = haplo.count$Count.Cases / Ncases *100, .after = "Count.Cases")
  haplo.count <- haplo.count %>% add_column(Freq.Controls = haplo.count$Count.Controls /Ncontrols *100, .after = "Count.Controls")
  
  # Create vector and initialize
  haplo.count$OR <- NA; haplo.count$Chi2 <- NA; As2control <- settings$allele2control %>% unlist()
  haplo.count$RefCases <- NA; haplo.count$RefControls <- NA
  
  # For each haplotype, compute Chi SQ based on reference 
  for (i in 1:nrow(haplo.count)){
    
    # Get ids of reference and filter out based on the allele to control / exclude
    HLA.ref <- HLA.df %>% filter(sample.id %notin% subj.ids.idx[hap1.code == i | hap2.code == i])
    for(A in As2control){
      locus <- A %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
      locusIds <- paste0(rep(locus,2), c(".1", ".2"))
      allele.id <-  A %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)
      HLA.ref <- HLA.ref %>% filter_at(locusIds, all_vars(.!=allele.id))
    }

    # Get ref number of cases and controls, append
    RefCases <- HLA.ref$pheno %>% table() %>% .["1"]; RefControls <- HLA.ref$pheno %>% table() %>% .["0"]
    haplo.count$RefCases[i] <- RefCases; haplo.count$RefControls[i] <- RefControls
    
    # Create contingency table
    cont.table <- matrix(c(haplo.count[i,"Count.Cases"], RefCases, haplo.count[i,"Count.Controls"], RefControls), ncol = 2)
    
    # Compute OR and Chi2
    chi2.res <- chisq.test(cont.table)
    OR.res <- (cont.table[1]*cont.table[4])/(cont.table[2]*cont.table[3])
    
    # Append to dataframe 
    haplo.count$Chi2[i] <- chi2.res$p.value
    haplo.count$OR[i] <- OR.res
  }
  
  # Return count of and stats of stats 
  return(haplo.count)
  
}

########### HAPLOTYPE ANALYSIS  ###########  

# Merge DRB3-4-5 into a single DRB locus
HLA.df <- merge_DRB(settings = settings, HLA.df)

# Check that DRB3-4-5 were merged and compute mean probability of loci for filtering
if (c("DRB.1", "DRB.2") %in% colnames(HLA.df) %>% any()){
  
  # Compute mean probabilities
  probs.df$prob.DRB <- apply(probs.df[c(grepl(colnames(probs.df), pattern = paste(c("DRB3", "DRB4", "DRB5"), collapse = "|")))], 
                             MARGIN = 1, function(x) mean(x))
  
}

# Count haplotypes with Expectation Maximization and run Chi2 test
haplo.count <- haplo_EM(settings, HLA.df, probs.df)

# Write 
write.xlsx(x = haplo.count, file = paste0(settings$Output$Haplotype, "HLA_EM.xlsx"), col.names = TRUE, row.names = FALSE)

