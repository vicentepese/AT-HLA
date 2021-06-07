## ---------------------------
##
## Name: HLA_glm_iter.R
##
## Desciption: Given a dataset of HLA calls, it fits GLM controlling
##                  for PCs and iterativelly controls for the most 
##                  significant allele.
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

########## IMPORT ########## 

# Import settings
settings <- jsonlite::read_json("settings.json")
options(stringsAsFactors = F)

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
as2control <- settings$allele2control

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


# Exclude allele
allele2exclude <- settings$allele2exclude %>% unlist()
a2exclude <- settings$allele2exclude %>% unlist();
if (!a2exclude %>% is_empty()){
  for (allele in a2exclude){
    
    # Parse locus and allele
    locus <- allele %>% strsplit("\\*") %>% unlist() %>% head(n=1)
    A <- allele %>% strsplit("\\*") %>% unlist() %>% tail(n=1)
    
    # Filter HLA calls 
    HLA.df <- HLA.df[which(HLA.df[,paste0(locus,".1")] != A & HLA.df[,paste0(locus,".2")] != A),]
    
  }
}

# Delete files to allow output to be written
file.names <- list.files(settings$Output$GLM, full.names = TRUE)
file.names <- file.names[grepl(file.names, pattern = "iter")]
file.remove(file.names)

########## ONE HOT ENCODING FUNCTIONS ########## 

# OHE for allele
alleleFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id~make_HLA, fun.aggregate = length)
  dcast_HLA$sample.id = as.factor(dcast_HLA$sample.id)
  return(setDF(dcast_HLA))
}

# OHE for carrier
carrierFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id~make_HLA, fun.aggregate = length)
  dcast_vals <- dcast_HLA[,!c("sample.id")]; dcast_vals[dcast_vals > 1] <- 1
  dcast_HLA <- cbind(data.frame(sample.id =dcast_HLA$sample.id), dcast_vals)
  dcast_HLA$sample.id = as.factor(dcast_HLA$sample.id)
  return(setDF(dcast_HLA))
}

########## ALLELE AND CARRIER FREQUENCY ########## 

computeACFREQ = function(data, locus, Dx){

  # Parse alleles
  A1 <- locus %>% paste0('.1'); A2 <- locus %>% paste0('.2')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
  carriers <- data[,c(A1, A2)]
  carriers.levels <- list(data[, A1], 
                          data[, A2]) %>% unlist() %>% levels()
  
  # Compute Carrier counts and frequencies
  carriers.unique <- apply(carriers, 1, function(x) unique(x)) %>% unlist() %>% as.factor()
  carriers.count <- table (carriers.unique); carriers.count[c(carriers.levels %>% setdiff(carriers.count %>% names()))] <- 0
  carriers.freq <- carriers.count /nrow(data) * 100
  carrier.df <- data.frame(allele = carriers.count %>% names(), 
                           carrierCount = carriers.count %>% as.vector(),
                           carrierFreq = carriers.freq %>% as.vector(),
                           carrierTotal = nrow(data))
  
  # Heterozigous, homozigous, and absence count
  A0 <- c(); A1 <- c(); A2 <- c();
  for (A in levels(as.factor(alleles))){
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  
  # Replace NAs with 0s for counts
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(as.factor(alleles)), A0 = A0, A1 = A1, A2 = A2)
  
  # Merge 
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele')
  
  # Rename columns and return based on diagnosis
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Case',4), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Control',4), sep = ''))
            return(ACFREQ.df)
          }
  )
}

########## CONTROL ALLELES ########## 

controlAllele = function(as2control, HLA.df){
  
  # Get unique loci 
  lociControl = lapply(as2control, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist() %>% unique() 
  
  # For each allele 
  alleleControl.df = data.frame(sample.id = HLA.df$sample.id)
  for (A in as2control){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    allele2control = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    allele1 <- paste(locus, '.1', sep = '')
    allele2 <- paste(locus, '.2', sep = '')
    
    # Get subjects
    alleleControl.df[A] <- as.logical(c(HLA.df[,c(allele1)] %>% as.character() == allele2control) + 
                                        c(HLA.df[,c(allele2)] %>% as.character() == allele2control)) %>% as.integer()
    
  }
  
  # If dataframe is empty, create 
  if(alleleControl.df %>% is_empty()){
    alleleControl.df <- data.frame(sample.id = HLA.df$sample.id)
  }
  
  return(alleleControl.df)
}


############ REGRESSION MODEL ############

runLogisticRegression = function(locus, OHE.carrierFreq.data, covars.df, as2control = NULL){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  alleles.freq <- colnames(OHE.carrierFreq.data)[-c(1,(ncol(OHE.carrierFreq.data)-length(as2control)):ncol(OHE.carrierFreq.data))]
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, covars.df[,c("sample.id", "PC1", "PC2","PC3")], by = 'sample.id')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = locus)]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  
  # For each alleles
  for (allele in alleles.freq){
    
    # Create string corresponding to alleles to control
    if (!is_empty(as2control)){
      control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.alleles <- ''
    }
    
    # Create GLM formula as string
    glm.formula <- paste('pheno ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    
    # Fit GLM model
    Acarrier.model <- glm(data = OHE.carrierFreq.data, 
                          formula = as.formula(glm.formula),
                          family = 'binomial', maxit = 100) %>% summary()
    
    # Format output
    Acarrier.model.df <- rbind(Acarrier.model.df, c(Acarrier.model$coefficients[2,1], 
                                                    Acarrier.model$coefficients[,dim(Acarrier.model$coefficients)[2]]))
  }
  
  # Change column names of output
  colnames(Acarrier.model.df) <- c('allele.COEF.CARRIER', c('Incercept', 'allele', Acarrier.model$coefficients[-c(1,2),] %>% row.names()) %>%
                                     paste('.CARRIER.pval', sep = ''))
  
  # Merge to include Alleles
  Acarrier.model.df <- data.frame(allele=alleles.freq, Acarrier.model.df)

  # Return
  return(Acarrier.model.df)
  
}

fitGLM = function(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control = NULL){
  
  # Control for allele
  controlAllele.df <- controlAllele(as2control, HLA.df)
  
  # Subset locus
  allele1 <- paste(locus, '.1',sep = '')
  allele2 <- paste(locus,'.2', sep =  '')
  data.locus <- HLA.df[,c('sample.id', allele1, allele2)] %>%
    filter(sample.id %in% c(data.cases$sample.id, data.controls$sample.id))
  
  # Compute allele frequencies and coutns, and carrier frequencies and counts
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele %notin% c('',"00:00"))
  
  # Parse one hot encoding and merge
  OHE.carrierFreq.data <- carrierFreqOHE(data.locus); 
  if ('V1' %in% colnames(OHE.carrierFreq.data)){OHE.carrierFreq.data$V1 <- NULL}
  OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), 
                                covars.df[c('pheno', 'sample.id')], by.x ='sample.id', by.y = 'sample.id') %>% 
    merge(controlAllele.df, by = 'sample.id')
  
  # Remove subjects thar are not controls or cases
  OHE.carrierFreq.data <- OHE.carrierFreq.data %>% filter(pheno != -9)
  
  # Remove useless alleles
  OHE.carrierFreq.data <- OHE.carrierFreq.data[,which(names(OHE.carrierFreq.data) %notin% c("", "00:00"))]
  
  # Run logistic regression model 
  glm.data <- runLogisticRegression(locus, OHE.carrierFreq.data, covars.df, as2control)
  
  # Create dataframes 
  HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                              ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                              by = 'allele')
  
  return(HLA.GLM_carriers.df)
}

############### HLA ANALYSIS ###########


# Get cases and controls
cases.ids <- covars.df$sample.id[which(covars.df$pheno ==1)]
controls.ids <- covars.df$sample.id[which(covars.df$pheno ==0)]
data.cases <- HLA.df %>% filter(sample.id %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id %in% controls.ids)

# Initialize while lopp
pval <- 0; signAlleles <- list(); 

# HLA Loci
loci <- colnames(HLA.df)[grepl(colnames(HLA.df), pattern = "\\.1")]
loci <- loci %>% lapply(function(x) strsplit(x, split = "\\.") %>% unlist() %>% head(n=1)) %>% unlist() 

# While signifiant alleles
idx <- 1; 
while (pval < 0.05){
  
  # Verbose:
  print(paste("Iteration", as.character(idx)))
  
  # Fit GLM to each allele of each locus
  HLA.GLM_carriers.list <- list(); pvalTotal <- c()
  for (locus in loci){
    
    # Verbose
    print(paste("Current locus: HLA-", locus, sep = ""))
    
    # Filter out subjects with imputation probability threshold
    probs.df_filt <- probs.df %>% filter(get(paste0("prob.", locus)) > prob_thr)
    data.cases.filt <- data.cases %>% filter(sample.id %in% probs.df_filt$sample.id)
    data.controls.filt <- data.controls %>% filter(sample.id %in% probs.df_filt$sample.id)
    
    # Fit GLM
    HLA.GLM_carriers.df <- fitGLM(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control)
    
    # Filter out low frequencies
    HLA.GLM_carriers.df_filt <- HLA.GLM_carriers.df %>% filter(carrierFreqCase > freq_thr | carrierFreqControl > freq_thr)
    
    # Apply pvalue correction
    HLA.GLM_carriers.df_filt <- add_column(HLA.GLM_carriers.df_filt, 
                                           allele.CARRIER.pval.CORR = p.adjust(p = HLA.GLM_carriers.df_filt$allele.CARRIER.pval, method = "BY"), 
                                           .after = "allele.CARRIER.pval")
    
    # Rbind with non-corrected alleles 
    HLA.GLM_carriers.df_filt <- HLA.GLM_carriers.df_filt %>% rbind.fill(HLA.GLM_carriers.df %>% filter(allele %notin% HLA.GLM_carriers.df_filt$allele))
    
    # Skip iteration to avoid colinearity
    for (skipAllele in settings$skipAllele %>% unlist()){
      L2skip <- skipAllele %>% strsplit(split = "\\*") %>% unlist() %>% head(n=1)
      if (locus == L2skip){
        A2skip <- skipAllele %>% strsplit(split = "\\*") %>% unlist() %>% tail(n=1)
        HLA.GLM_carriers.df_filt <- HLA.GLM_carriers.df_filt[-which(HLA.GLM_carriers.df_filt$allele == A2skip),]
      }
    }
    
    # Append to list
    HLA.GLM_carriers.list[[locus]] <- HLA.GLM_carriers.df_filt 
    for (A in HLA.GLM_carriers.df_filt$allele){
      pvalA <- HLA.GLM_carriers.df_filt$allele.CARRIER.pval.CORR[which(HLA.GLM_carriers.df_filt$allele == A)]
      pvalTotal[paste(locus, '*', A, sep = '')] <- pvalA
    }
  }
  
  # Exlude NAs
  pvalTotal <- pvalTotal[!is.na(pvalTotal)]

  # Get minimum allele value 
  pval <- pvalTotal[which(pvalTotal == min(pvalTotal))][1]
  pvalMin <- pval %>% names(); pvalMinLocus <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[1]
  pvalMinAllele <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[2]
  
  # Add to outputs
  as2control <- c(as2control, pvalMin)
  preOut <- HLA.GLM_carriers.list[[ pvalMinLocus]] %>% filter(allele == pvalMinAllele)
  preOut$allele <- pvalMin  ; signAlleles[[idx]] <- preOut
  
  # Update index
  idx <- idx +1
  
}

# Format output
allele <- as.data.frame(signAlleles[[1]]); 
if (length(signAlleles)>1){
  for (idx in 2:length(signAlleles)){
    allele <- rbind.fill(allele, signAlleles[[idx]] %>% as.data.frame())
  }
}

# Write outcome
for (idx in 1:length(HLA.GLM_carriers.list)){
  HLA.GLM_carriers.df <- HLA.GLM_carriers.list[[idx]]; locus <- HLA.GLM_carriers.list %>% names() %>% .[idx]
  write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$Output$GLM, 'HLA_GLM_Carriers_iter','.xlsx', sep = ''), sheetName = locus, 
             col.names = TRUE, row.names = FALSE, append = TRUE)
}
write.xlsx(x = allele, file = paste(settings$Output$GLM, 'HLA_GLM_Carriers_iter','.xlsx', sep = ''), sheetName = 'Significant_alleles', 
           col.names = TRUE, row.names = FALSE, append = TRUE)



