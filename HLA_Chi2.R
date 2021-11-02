## ---------------------------
##
## Name: HLA_analysis.R
##
## Desciption: Given a dataset of HLA calls, computes a Chi2 
##    for each locus and allele between cases and controls.
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
library(jsonlite, warn.conflicts = F)
library(plyr, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse))
library(readr,warn.conflicts = F)
library(data.table, warn.conflicts = F)
library(xlsx,warn.conflicts = F)
library(utils,warn.conflicts = F)

########## ERROR HANDLING ########### 

settingsCheck = function(settings){
  
  # Check files existence
  if (!settings$file$HLA_Data %>% file.exists()){
    stop("HLA data file does not exist.")
  } else if (!settings$file$covars %>% file.exists()){
    stop("Covariates file does not exist.")
  } else if(!settings$file$probs %>% file.exists()){
    stop("Imputation probabilities file does not exist.")
  }
  
  # Check parameters 
  if (settings$prob_thr < 0 | settings$prob_thr > 1){
    stop("The probability threshold must be between 0 and 1.")
  }
  if (settings$freq_thr < 0 | settings$prob_thr >1){
    stop("The frequency threshold must be between 0 and 1.")
  }
  
  # Throw warning 
  if (settings$prob_thr == 1){
    warning("Probability threshold is equal to 1. All cases will be included.")
  }
  
}

phenoCheck = function(covars.df){

  # Get phenotype
  pheno_labels <- covars.df$pheno %>% unique()

  # Check that phenotype is between 0 and 1, or 1 and 2 
  if (any(pheno_labels %notin% c(1,2))){  
    covars.df$pheno <- covars.df$pheno - 1
  } else if (any(pheno_labels %notin% c(0,1))){ 
    stop("Phenotype labels are not 1 and 2, or 0 and 1.")
  }
  
  # Get number of cases and controls
  pheno_count <- covars.df$pheno  %>% table()
  
  # Check that there are cases and controls
  if (pheno_count['0'] == 0){
    stop("No controls in the covariates file.")
  } else if (pheno['1'] == 0){
    stop("No cases in the covariates file.")
  }
  
  return(covars.df)
}

matchedControlsCheck = function(matched_controls){

  # Check that the list is not empty
  if (matched_controls %>% is_empty){
    stop ("No matched controls were found for the given cases.")
  }

}

alleleCheck = function(HLA.df, locus, A){

  # Get alleles
  alleles <- unique(c(HLA.df[,c(paste0(locus,".1"))], HLA.df[,c(paste0(locus,".2"))]))
  
  # Check that allele existst
  if (A %notin% alleles){
    stop("Allele to exclude not present in the data.")
  }

}

ethnicityCheck = function(settings, ethnicity.df){
  
  # Check that ethnicities are in the data provided 
  if (any(settings$ethnicity %notin% ethnicity.df$Population)){
    stop("One of the provided ethnicities is not present in the data.")
  }
  
}

########### INITIALIZATION ########### 

# Import settings
settings <- jsonlite::read_json("settings.json")
options(stringsAsFactors = F)

# Check settings
settingsCheck(settings)

# Verbose 
if (settings$verbose) print("Loading data, covariates, and imputation probabilities.")

# Create command
`%notin%` <- Negate(`%in%`)

# Import HLA calls, covariates, and matched controls
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Check pheno
covars.df <- phenoCheck(covars.df)

# If list of matched controls provided, filter
if (!settings$file$matched_controls %>% is_empty()){
  
  # Verbose
  if (settings$verbose) print("Parsing matched controls.")
  
  # Get cases ids, and cases from data (may not be the same, e.g. sub-dataset of only tumors)
  cases.ids <- covars.df %>% filter(pheno == 1) %>% select(sample.id) %>% unlist()
  HLA.cases.ids <- HLA.df %>% filter(sample.id %in% cases.ids) %>% .["sample.id"] %>% unlist()
  
  # Get list of matched controls and merge 
  match_cntrls.df <- read.csv(settings$file$matched_controls)
  matched_controls <- match_cntrls.df %>% 
    filter(sample.id_case %in% HLA.cases.ids) %>% .[,2:ncol(match_cntrls.df)] %>% flatten() %>% unlist() %>% unique()
  ids <- c(HLA.cases.ids, matched_controls)
  
  # Check matched controls 
  matchedControlsCheck(matched_controls)
  
  # Filter the HLA data, covariates and probabilities
  HLA.df <- HLA.df %>% filter(sample.id %in% ids)
  covars.df <- covars.df %>% filter(sample.id %in% ids)
  probs.df <- probs.df %>% filter(sample.id %in% ids)
}

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100
allele2exclude <- settings$allele2exclude %>% unlist()

# Filter out the alleles to exclude 
if (!allele2exclude %>% is_empty()){
  
  # Verbose 
  if (settings$verbose) print("Excluded alleles: ")
  
  for (allele in allele2exclude){

    # Parse locus and allele
    locus <- allele %>% strsplit("\\*") %>% unlist() %>% head(n=1)
    A <- allele %>% strsplit("\\*") %>% unlist() %>% tail(n=1)
    
    # Check allele
    alleleCheck(HLA.df, locus, A)
    
    # Filter HLA calls 
    HLA.df <- HLA.df[which(HLA.df[,paste0(locus,".1")] != A & HLA.df[,paste0(locus,".2")] != A),]
    
    # Verbose
    if (settings$verbose) print(allele2exclude)
    
  }
}

# Parse HLA calls based on ethnicity, if provided
if (!settings$ethnicity %>% is_empty()){
  
  # Parse IDs in ethnicity/ies
  ethnicity.df <- read.csv(settings$file$ethnicity)
  
  # Check ethnicity
  ethnicityCheck(settings, ethnicity.df)
  
  # Filter
  ethnicity.df.filt <- ethnicity.df %>% filter(Population %in% settings$ethnicity %>% unlist())
  HLA.df <- HLA.df %>% 
    filter(sample.id %in% ethnicity.df.filt$sample.id)
  
  # Verbose
  if (settings$verbose) print(paste("Ethnicities included:", paste(settings$ethnicity, sep =" "), sep = " "))

}

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

# Verbose 
if (settings$verbose) print("Deleting previous files.")

# Delete files to allow output to be written
file.names <- list.files(settings$Output$Chi2, full.names = TRUE)
file.remove(file.names)

########### COMPUTE ALLELE/CARRIER COUNT/FREQUENCIES ########### 

computeACFREQ = function(data, locus, Dx){
  # Compute heterozigous, homozigous, or absence count. 
  # Compute allele frequency, count and total.
  # Compute carrier frequency, count and total.
  
  # Allele frequency 
  A1 <- locus %>% paste0('.1'); A2 <- locus %>% paste0('.2')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  alleles.df <- data.frame(allele = alleles.count %>% names(), 
                           alleleCount = alleles.count %>% as.vector(),
                           alleleFreq = alleles.freq %>% as.vector(), 
                           alleleTotal = nrow(data)*2)
  
  # Carrier frequency
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
  
  # Heterozigous, homozigous, and absence count
  A0 <- c(); A1 <- c(); A2 <- c();
  for (A in levels(as.factor(alleles))){
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(as.factor(alleles)), A0 = A0, A1 = A1, A2 = A2)
  
  # Merge 
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele') %>% merge(alleles.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',7), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',7), sep = ''))
            return(ACFREQ.df)
          }
  )
}


########### COMPUTE CHI SQUARE ########### 

computeChi2 = function(ACFREQ.df){
  
  # Initialize 
  chi2.pval.alleles <- c() ; chi2.pval.carriers <- c();
  chi2.OR.alleles <- c(); chi2.OR.carriers <- c();
  fishers.pval.alleles <- c(); fishers.pval.carriers <- c()
  fishers.OR.alleles <- c(); fishers.OR.carriers <- c()
  fishers.UI.alleles <- c(); fishers.UI.carriers <- c()
  fishers.LI.alleles <- c(); fishers.LI.carriers <- c()
  OR.carrier <- c(); OR.allele <- c()
  
  # For each allele in the locus 
  for (A in ACFREQ.df$allele){
    ## Allele frequency 
    # Create contingency table  and run chi2 test
    allele.data <- ACFREQ.df %>% filter(allele == A)
    cont.table.allele <- matrix(c(allele.data$alleleCountCase,
                                  allele.data$alleleCountControl,
                                  allele.data$alleleTotalCase - allele.data$alleleCountCase,
                                  allele.data$alleleTotalControl - allele.data$alleleCountControl),
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
    # OR 
    OR.allele <- c(OR.allele, (cont.table.allele[1]*cont.table.allele[4])/(cont.table.allele[2]*cont.table.allele[3]))
    
    
    ## Carrier Frequency
    # Create contingency table  and run chi2 test
    cont.table.carrier <- matrix(c(allele.data$carrierCountCase,
                                   allele.data$carrierCountControl,
                                   allele.data$carrierTotalCase - allele.data$carrierCountCase,
                                   allele.data$carrierTotalControl - allele.data$carrierCountControl),
                                 nrow = 2)
    # Tests
    chi2.carrier.res <- chisq.test(x = cont.table.carrier)
    fishers.carrier.res <- fisher.test(cont.table.carrier)
    # Chi2
    chi2.pval.carriers <- c(chi2.pval.carriers, chi2.carrier.res$p.value) ; 
    chi2.OR.carriers <- c(chi2.OR.carriers, chi2.carrier.res$statistic %>% unname())
    # Fishers' exact test 
    fishers.pval.carriers <- c(fishers.pval.carriers, fishers.carrier.res$p.value);
    fishers.OR.carriers <- c(fishers.OR.carriers, fishers.carrier.res$estimate %>% unname())
    fishers.LI.carriers <- c(fishers.LI.carriers, fishers.carrier.res$conf.int[1])
    fishers.UI.carriers <- c(fishers.UI.carriers, fishers.carrier.res$conf.int[2])
    # OR 
    OR.carrier <- c(OR.carrier, (cont.table.carrier[1]*cont.table.carrier[4])/(cont.table.carrier[2]*cont.table.carrier[3]))
    
  }
  
  # Create dataframe 
  stats.data <- data.frame(allele = ACFREQ.df$allele, FishersCarrierPVAL = fishers.pval.carriers, 
                           FishersCarrierOR = fishers.OR.carriers, FishersCarrierLI = fishers.LI.carriers,
                           FishersCarrierUI = fishers.UI.carriers,
                           ChiCarrierPVAL = chi2.pval.carriers, ChiCarrierEST = chi2.OR.carriers, ORCarrier = OR.carrier,
                           FishersAllelePVAL = fishers.pval.alleles, FishersAlleleOR = fishers.OR.alleles, 
                           FishersAlleleLI = fishers.LI.alleles, FishersAlleleUI = fishers.UI.alleles,
                           ChiAllelePVAL = chi2.pval.alleles, ChiAlleleEST = chi2.OR.alleles,
                           ORAllele = OR.allele)
  
  # Remove NAs
  stats.data[is.na(stats.data)] <- 0
  
  # Return 
  return(stats.data)
  
}


########### MAIN LOOP ########### 

# Get cases and controls and separate datasets
cases.ids <- covars.df$sample.id[which(covars.df$pheno == 1)]
controls.ids <- covars.df$sample.id[which(covars.df$pheno == 0)]
data.cases <- HLA.df %>% filter(sample.id %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id %in% controls.ids)

# Initialize multiple-test correction 
pvals <- c(); l1group <- c(); l2group <- c(); l2locus <- c()

# For each locus 
loci <- colnames(HLA.df)[grepl(colnames(HLA.df), pattern = "\\.1")]
loci <- loci %>% lapply(function(x) strsplit(x, split = "\\.") %>% unlist() %>% head(n=1)) %>% unlist() 
idx <- 1

for (locus in loci){
  
  # Verbose
  print(paste0("Conducting analysis on locus: ", locus))
  
  # Filter out subjects with imputation probability threshold
  probs.df_filt <- probs.df %>% filter(get(paste0("prob.", locus)) > prob_thr)
  data.cases.filt <- data.cases %>% filter(sample.id %in% probs.df_filt$sample.id)
  data.controls.filt <- data.controls %>% filter(sample.id %in% probs.df_filt$sample.id)
  
  # Compute allele and carrier counts and frequencies
  ACFREQ.cases <- computeACFREQ(data.cases.filt, locus, 'case');
  totalCases <-unique(ACFREQ.cases$alleleTotalCase); carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls.filt, locus, 'control');
  totalControls <-unique(ACFREQ.controls$alleleTotalControl); carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean 
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df)); ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df));
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele %notin% c('','00:00'))
  
  # Compute Chi2
  stats.data <- computeChi2(ACFREQ.df)
  stats.data[stats.data==Inf] = 'Inf'
  
  # Create dataframes 
  HLA.alleles.df <- merge(stats.data[,c(1,which(grepl('Allele', colnames(stats.data))))],
                          ACFREQ.df[,c(which(grepl(paste(c('A0','A1','A2', 'allele'), collapse = '|'), colnames(ACFREQ.df))))],
                          by = 'allele')
  HLA.carriers.df <-  merge(stats.data[,c(1,which(grepl('Carrier', colnames(stats.data))))],
                            ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                            by = 'allele')
  
  # Filter out low frequencies 
  HLA.alleles.df_filt <- HLA.alleles.df %>% filter(alleleFreqCase > freq_thr | alleleFreqControl > freq_thr)
  HLA.carriers.df_filt <- HLA.carriers.df %>% filter(carrierFreqCase > freq_thr | carrierFreqControl > freq_thr)
  
  # Apply p-value correction
  HLA.alleles.df_filt <- add_column(HLA.alleles.df_filt,
                               FishersAllelePVAL_CORR = p.adjust(p = HLA.alleles.df_filt$FishersAllelePVAL, method = "BY"),
                               .after = "FishersAllelePVAL")
  HLA.alleles.df_filt <- add_column(HLA.alleles.df_filt,
                               ChiAllelePVAL_CORR = p.adjust(p = HLA.alleles.df_filt$ChiAllelePVAL, method = "BY"),
                               .after = "ChiAllelePVAL")
  HLA.carriers.df_filt <- add_column(HLA.carriers.df_filt,
                               FishersCarrierPVAL_CORR = p.adjust(p = HLA.carriers.df_filt$FishersCarrierPVAL, method = "BY"),
                               .after = "FishersCarrierPVAL")
  HLA.carriers.df_filt <- add_column(HLA.carriers.df_filt,
                                ChiCarrierPVAL_CORR = p.adjust(p = HLA.carriers.df_filt$ChiCarrierPVAL, method = "BY"),
                                .after = "ChiCarrierPVAL")
  
  # Rbind with non-corrected alleles 
  HLA.alleles.df_filt <- HLA.alleles.df_filt %>% rbind.fill(HLA.alleles.df %>% filter(allele %notin% HLA.alleles.df_filt$allele))
  HLA.carriers.df_filt <- HLA.carriers.df_filt %>% rbind.fill(HLA.carriers.df %>% filter(allele %notin% HLA.carriers.df_filt$allele))
  
  # Write 
  write.xlsx(x = HLA.alleles.df_filt, file = paste0(settings$Output$Chi2, 'HLA_AnalysisAlleles','.xlsx', sep = ''), sheetName = locus,
             col.names = TRUE, row.names = FALSE, append = TRUE)
  write.xlsx(x = HLA.carriers.df_filt, file = paste0(settings$Output$Chi2, 'HLA_AnalysisCarriers','.xlsx', sep = ''), sheetName = locus,
             col.names = TRUE, row.names = FALSE, append = TRUE)
  
}

# Verbose
if (settings$verbose) print(paste("Outputs saved in:", settings$Output$Chi2, sep=" "))

