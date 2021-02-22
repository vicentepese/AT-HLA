## ---------------------------
##
## Name: HLA_AA_glm.R
##
## Desciption: Given a dataset of HLA calls, it tests the significance
##              of aminoacid-position pairs based on each allele by 
##              fitting a GLM for each pair controlling for PCs.
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
library(TreeBH)
library(zeallot)
library(epitools)

########## INITIALIZATION #########

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

# Read AAs to control
AAs2control <- settings$AA2control

########## OHE FUNCTIONS ##########

AA2OHE = function(settings, AAs, pos, L, AA_locus, HLA.df){
  
  # Get locus ID
  c(A1, A2) %<-% c(paste0(L,'.1'), paste0(L,'.2'))
  
  # Initialize AAs in OHE format 
  OHE.AA <- data.frame(sample.id = HLA.df$sample.id)
  
  # For each AA
  for (AA in AAs){
    
    # Get alleles with AAs
    allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
    allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                           .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
    
    # Get data with such alleles
    OHE.AA[paste(L,pos,AA, sep = "_")] <- apply(HLA.df[,c(A1,A2)], 2, function(x,allelesAA) as.integer(x %in% allelesAA), allelesAA) %>%
                           apply(1,function(x) as.integer(x[1] | x[2]))
  }
 
  # Return 
  return(OHE.AA)
   
}

controlAAs = function(AAs2control, HLA.df){
  
  # For each AA
  
}

########## COUNT AMINOACIDS ############

count_AA = function(OHE.AA, locus, pos, AAs){
  
  # Parse cases and controls 
  OHE.AA.cases <- OHE.AA %>% filter(pheno == 1)
  OHE.AA.controls <- OHE.AA %>% filter(pheno == 0)
  
  # Parse AA 
  AAs2control <- settings$AA2control %>% unlist()
  AAs.IDs <- colnames(OHE.AA)[-c(1,(ncol(OHE.AA)-length(AAs2control)):ncol(OHE.AA))]
  
  # Count cases for each allele 
  Ncases <- OHE.AA.cases[,AAs.IDs] %>% apply(MARGIN=2, function (x) x %>% table() %>% .['1'])
  Ncontrols <- OHE.AA.controls[,AAs.IDs] %>% apply(MARGIN = 2, function (x) x %>% table() %>% .['1'])
  
  # Calculate frequencies 
  FreqCases <- Ncases/nrow(OHE.AA.cases) *100
  FreqControls <- Ncontrols/nrow(OHE.AA.controls) * 100
  
  # Create dataframe 
  AA.count <- data.frame(locus = rep(locus, length(AAs)), pos = rep(pos, length(AAs)), AA = AAs, 
                         Ncases = Ncases, FreqCases = FreqCases,
                         Ncontrols = Ncontrols, FreqControls = FreqControls)
  
  # Return 
  return(AA.count)
  
}

########### PARSE ALLELES ##########

parseAlleles = function(AA.model.df, AA_locus){
  
  # For each Amino Acid, get the alleles where it is present 
  alleles <- list()
  for (i in 1:nrow(AA.model.df)){
    pos <- AA.model.df[i,"pos"]; AA <- AA.model.df[i,"AA"]
    alleles[[i]] <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
  }
  
  # Append to dataframe 
  AA.model.df$alleles <- alleles %>% lapply(function(x) x %>% unname() %>% paste(collapse = ", ")) %>% unlist()
  
  # Return 
  return(AA.model.df)

  
}

############ REGRESSION MODEL ############

run_GLM = function(OHE.AA, covars.df, L, pos){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  AAs2control <- settings$AA2control %>% unlist()
  AAs.IDs <- colnames(OHE.AA)[-c(1,(ncol(OHE.AA)-length(AAs2control)):ncol(OHE.AA))]
  OHE.AA <- merge(OHE.AA, covars.df[,c("sample.id", "PC1", "PC2", "PC3")], by = 'sample.id')
  
  # Remove alleles for control
  OHE.AA[AAs2control] <- NULL
  AAs.IDs <- AAs.IDs[!AAs.IDs %in% AAs2control]
  
  # Run logistic regression on carrier frequency 
  AA.model.df <- data.frame(); AAs <- c()
  for (AA in AAs.IDs){
    
    # Ensure presence of allele / skip aliased 
    if (OHE.AA[, AA] %>% unique() %>% length() == 1){
      next
    }
    
    # Append AA
    AAs <- c(AAs, AA %>% strsplit("_") %>% unlist() %>% tail(n=1))
    
    # Control for AA (if present, else empty)
    if (length(AAs2control) > 0){
      control.AAs <- paste(' ', AAs2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.AAs <- ''
    }
    
    # Fit GLM and bind row-wise to intial dataframe 
    glm.formula <- paste('pheno ~ `',AA, '` + PC1 + PC2 + PC3', control.AAs, sep = '')
    AA.model <- glm(data = OHE.AA, 
                          formula = as.formula(glm.formula),
                          family = 'binomial', maxit = 100) %>% summary()
    AA.model.df <- rbind(AA.model.df, c(AA.model$coefficients[2,1], 
                                        AA.model$coefficients[,dim(AA.model$coefficients)[2]]))
    
  }
  
  # Change column names 
  colnames(AA.model.df) <- c('AA.COEF', c('Incercept', 'AA', AA.model$coefficients[-c(1,2),] %>% row.names()) %>%
                                     paste('.pval', sep = ''))
  
  # Create data frame 
  AA.model.df <- data.frame(locus = rep(L, nrow(AA.model.df)), pos = rep(pos, nrow(AA.model.df)), AA = AAs,
                            AA.model.df)

  # Return
  return(AA.model.df)
  
}

########## AMINO ACID ANALYSIS ##########

# Get loci, number of cases and controls
loci <- unique(AA_alignment$locus)
c(Ncases, Ncontrols) %<-% c(HLA.df$pheno %>% table %>%.['1'], HLA.df$pheno %>% table %>%.['0'])

# Initialize loop
AA.df <- data.frame()

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
      
      # Convert AA carriers into OHE
      OHE.AA <- AA2OHE(settings, AAs, pos, L, AA_locus, HLA.df)
      
      # Append phenotypes 
      OHE.AA <- merge(OHE.AA, covars.df[,c("sample.id", "pheno")], by = "sample.id")
      
      # Filter out subjects with low imputation probability
      probs.df_filt <- probs.df %>% filter(get(paste0("prob.", L)) > prob_thr)
      OHE.AA <- OHE.AA %>% filter(sample.id %in% probs.df_filt$sample.id)
      
      # Get OHE of control AAs 
      # OHE.controlAA <- AA2OHE(settings, AAs)
      
      # Fit GLM 
      AA.model.df <- run_GLM(OHE.AA, covars.df, L, pos)
      
      # Count
      AA.count <- count_AA(OHE.AA, L, pos, AAs)
      
      # Merge
      AA.model.df <- merge(AA.model.df, AA.count, by = c("locus", "pos", "AA"))
      
      # Parse alleles
      AA.model.df <- parseAlleles(AA.model.df, AA_locus)
      
      # Bind row-wise to dataframe 
      AA.df <- rbind(AA.df, AA.model.df)
    }
  }
}


# Adjust p-value through FDR correction, append to results
pval.corr <- p.adjust(p = AA.df$AA.pval, method = "BY")
AA.df <- add_column(AA.df, pval.corr, .after = "AA.pval")

# Save dataframe 
write.csv(AA.df, file = paste0(settings$Output$AA_GLM, "AA_GLM.csv"), row.names = FALSE)


