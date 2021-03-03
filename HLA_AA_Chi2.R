## ---------------------------
##
## Name: HLA_AA_Chi2.R
##
## Desciption: Given a dataset of HLA calls, it tests the significance
##              of aminoacid-position pairs based on each allele.
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
library(zeallot)
library(epitools)

########### INITIALIZATION ########### 

# Set working directory
setwd("~/Documents/HLA_association_pipeline")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls, covariates and probabilities
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Correct pheno
if (2 %in% covars.df$pheno %>% table() %>% names()){
  covars.df$pheno <- covars.df$pheno - 1
}

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

# Append pheno to HLA calls 
HLA.df <- merge(HLA.df, covars.df[,c("sample.id", "pheno")], by = "sample.id")

# Delete files to allow output to be written
file.names <- list.files(settings$Output$AA, full.names = TRUE)
file.remove(file.names)

# Import amino acid alignment 
AA_alignment <- read.table(settings$file$AA_alignment, header = TRUE, sep = ',')

# Convert from full protein to mature protein indexing
AA_alignment[,3] <- AA_alignment %>% apply(MARGIN = 1, function(x) x[3] %>% substr(start = 30, stop = nchar(x[3])))

########### EXTRA FUNCTIONS ########### 

carrierCount = function(settings, data){
  
  # Get loci, number of cases, number of controls
  loci <- colnames(data)[which(grepl('HLA', colnames(data)))] %>%
    lapply(function(x) x %>% strsplit('_') %>% unlist() %>% .[1] %>% strsplit('HLA') %>% unlist() %>% .[2]) %>% unlist() %>% unique()
  c(Ncases, Ncontrols) %<-% c(data$Dx %>% table() %>% .['1'], data$Dx %>% table() %>% .['0'])
  
  # Initialize loop 
  alleles <- c(); locusVec <- c(); cases.Count <- c(); controls.Count <- c()
  
  # For each locus, count carrier in cases and controls
  for (locus in loci){
    
    # Get alleles 
    c(A1, A2) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
    data.locus <- data[, c('GWASID', A1, A2, 'Dx')]
    
    # Get carriers and count 
    carriers.locus <- c(data.locus[,A1], data.locus[,A2]) %>% unique()
    carriers.locus <- carriers.locus[which(carriers.locus != "")]
    for(carrier in carriers.locus){
      data.carrier <- data.locus %>% filter(get(A1) == carrier | get(A2) == carrier)
      c(C.cases, C.controls ) %<-% c(table(data.carrier$Dx)['1'], table(data.carrier$Dx)['0']);
      if(C.cases %>% is.na()){C.cases <- 0}; if (C.controls %>% is.na()) {C.controls <- 0}
      c(F.cases, F.controls) %<-% c(C.cases/Ncases, C.controls/Ncontrols)
      
      # If frequency is smaller than a threshold, next. Else add to vector 
      if (is.na(C.cases) | is.na(C.controls)){
        next 
      } else if (F.cases < settings$min_carrierFreq | F.controls < settings$min_carrierFreq){
        next
      } else{
        locusVec <- c(locusVec, locus)
        A <- paste( locus, '*', carrier, ':01', sep = ''); alleles <- c(alleles, A)
        cases.Count <- c(cases.Count, C.cases); controls.Count <- c(controls.Count, C.controls)
      }
    }
  }
  
  # Create data frame 
  counts.data <- data.frame(locus = locusVec, allele = alleles, Ncases = cases.Count, Ncontrols = controls.Count)
  
  # Return 
  return (counts.data)
  
}

########### AMINO ACID ANALYSIS ########### 

# Get loci, number of cases and controls
loci <- unique(AA_alignment$locus)
c(Ncases, Ncontrols) %<-% c(HLA.df$pheno %>% table %>%.['1'], HLA.df$pheno %>% table %>%.['0'])

# Initialize loop
AA.total <- c(); pval <- c(); OR <- c(); pos.total <- c(); locus <- c();
cases.Count <- c(); controls.Count <- c(); alleles <- c();
cases.Freq <- c(); controls.Freq <- c();

# For each locus 
for (L in loci){
  
  # Verbose
  print(paste0("Current locus: HLA-", L))
  
  # Get locus ID
  c(A1, A2) %<-% c(paste0(L,'.1'), paste0(L,'.2'))
  
  # Get alignment subset, and counts subset
  AA_locus <- AA_alignment %>% filter(locus == L)
  
  # Get max sequence length 
  maxLen <- AA_locus$sequence %>% lapply(nchar) %>% unlist() %>% max()
  
  # Get unique AAs 
  for (pos in 1:maxLen){
    
    # Get unique AAs
    AAs <- AA_locus$sequence %>% lapply(function(x, pos) substr(x, pos, pos), pos) %>% unlist() %>% unique()
    
    # If more than one, count 
    if (length(AAs) >1){
      
      for (AA in AAs){
        
        # Get alleles with AAs
        allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
        allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                               .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
        
        # Get data with such alleles
        data.AA <- HLA.df %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
        
        # Filter out subjects with low imputation probability
        probs.df_filt <- probs.df %>% filter(get(paste0("prob.", L)) > prob_thr)
        data.AA <- data.AA %>% filter(sample.id %in% probs.df_filt$sample.id)
        
        # Filter out subjects for with low imputation probability in the total count
        HLA.df_filt <- HLA.df %>% filter(sample.id %in% probs.df_filt$sample.id)
        c(Ncases, Ncontrols) %<-% c(HLA.df_filt$pheno %>% table %>%.['1'], HLA.df_filt$pheno %>% table %>%.['0'])
        
        # Get counts and frequencies
        c(AA.cases, AA.controls) %<-% c(data.AA$pheno %>% table() %>% .['1'], data.AA$pheno %>% table() %>% .['0']);
        if(is.na(AA.cases)) {AA.cases <- 0}; if(is.na(AA.controls)){AA.controls <- 0}
        AA.casesFreq <- (AA.cases/Ncases)*100; AA.controlsFreq <- (AA.controls/Ncontrols)*100
        
        # Contingency table 
        cont.table <- matrix(c(AA.cases, AA.controls, 
                               Ncases - AA.cases, Ncontrols - AA.controls), nrow = 2)
        
        # Compute chiSq test
        chi.sq.res <- fisher.test(cont.table)
        
        # Append to vectors
        AA.total <- c(AA.total, AA); pval <- c(pval, chi.sq.res$p.value);
        OR <- c(OR, (cont.table[1]*cont.table[4])/(cont.table[2]*cont.table[3])); 
        pos.total <- c(pos.total, pos); locus <- c(locus, L);
        cases.Count <- c(cases.Count, AA.cases); controls.Count <- c(controls.Count, AA.controls)
        cases.Freq <- c(cases.Freq, AA.casesFreq); controls.Freq <- c(controls.Freq, AA.controlsFreq)
        alleles <- c(alleles, paste(allelesAA.OG, collapse = ', '))
      }
    }
  }
}

# Create dataframe
AA.analysis.results <- data.frame(locus= locus, AA = AA.total, pos = pos.total, Ncases = cases.Count, FreqCases = cases.Freq, 
                                  Ncontrol = controls.Count, FreqControls = controls.Freq,
                                  pval = pval, OR = OR, alleles = alleles)

# Adjust p-value through FDR correction, append to results
pval.corr <- p.adjust(p = AA.analysis.results$pval, method = "BY")
AA.analysis.results <- add_column(AA.analysis.results, pval.corr, .after = "pval")
  
# Save dataframe 
write.csv(AA.analysis.results, file = paste0(settings$Output$AA, "AA_Chi2.csv"), row.names = FALSE)


