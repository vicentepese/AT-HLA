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

########### INITIALIZATION ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls and covariates
if(settings$dataset == "Stanford"){
  HLA.df <- read.csv(settings$file$HLA_df)
  covars.df <- read.csv(settings$file$HLA_covars_filt)
  covars.df$pheno <- covars.df$pheno -1
  covars.df$sex <- covars.df$sex -1
} else if (settings$dataset == "UKB"){
  HLA.df <- read.csv(settings$file$HLA_df_UKB)
  covars.df <- read.csv(settings$file$HLA_covars_UKB)
  covars.df$pheno <- covars.df$pheno 
  
}



########### ONE HOT ENCODING FUNCTIONS ###############

# HLA Parse function
alleleFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id.file = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id.file, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id.file~make_HLA, fun.aggregate = length)
  dcast_HLA$sample.id.file = as.factor(dcast_HLA$sample.id.file)
  return(setDF(dcast_HLA))
}

# HLA Parse function
carrierFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id.file = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id.file, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id.file~make_HLA, fun.aggregate = length)
  dcast_vals <- dcast_HLA[,!c("sample.id.file")]; dcast_vals[dcast_vals > 1] <- 1
  dcast_HLA <- cbind(data.frame(sample.id.file =dcast_HLA$sample.id.file), dcast_vals)
  dcast_HLA$sample.id.file = as.factor(dcast_HLA$sample.id.file)
  return(setDF(dcast_HLA))
}

########### CONTROL ALLELES ##########

controlAllele = function(settings, locus, data.filt){
  
  # Alleles to control
  As2control = settings$controlAlleles
  
  # Get unique loci 
  lociControl = lapply(As2control, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist() %>% unique() 
  
  # For each allele 
  alleleControl.df = data.frame(sample.id.file = data.filt$sample.id.file)
  for (A in As2control){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    allele2control = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    allele1 <- paste(locus, '.1', sep = '')
    allele2 <- paste(locus, '.2', sep = '')
    
    # Get subjects
    alleleControl.df[A] <- as.logical(c(data.filt[,c(allele1)] %>% as.character() == allele2control) + 
                                        c(data.filt[,c(allele2)] %>% as.character() == allele2control)) %>% as.integer()
    
  }
  return(alleleControl.df)
}

################ COMPUTE ACFREQ ###############

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

############ REGRESSION MODEL ############


runLogisticRegression = function(locus, OHE.alleleFreq.data, OHE.carrierFreq.data, covars.df){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  as2control <- settings$controlAlleles %>% unlist()
  alleles.freq <- colnames(OHE.alleleFreq.data)[-c(1,(ncol(OHE.alleleFreq.data)-length(as2control)):ncol(OHE.alleleFreq.data))]
  OHE.alleleFreq.data <- merge(OHE.alleleFreq.data, covars.df[,-3], by = 'sample.id.file')
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, covars.df[,-3], by = 'sample.id.file')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = locus)]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.alleleFreq.data[allele.subset] <- NULL
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on allele frequency data
  Afreq.model.df <- data.frame()
  for (allele in alleles.freq){
    if (length(as2control) > 0){
      control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.alleles <- ''
    }
    glm.formula <- paste('pheno ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    Afreq.model <- glm(data = OHE.alleleFreq.data, 
                       formula = as.formula(glm.formula),
                       family = 'binomial', maxit = 100) %>% summary()
    Afreq.model.df <- rbind(Afreq.model.df, c(Afreq.model$coefficients[2,1], 
                                              Afreq.model$coefficients[,dim(Afreq.model$coefficients)[2]]))
  }
  colnames(Afreq.model.df) <- c('allele.COEF.ALLELE', 
                                c('Incercept', 'allele', Afreq.model$coefficients[-c(1,2),] %>% row.names()) %>% paste('.ALLELE.pval', sep = ''))
  Afreq.model.df <- cbind(data.frame(allele=alleles.freq, Afreq.model.df))
  
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  for (allele in alleles.freq){
    if (length(as2control) > 0){
      control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.alleles <- ''
    }
    glm.formula <- paste('pheno ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    Acarrier.model <- glm(data = OHE.carrierFreq.data, 
                          formula = as.formula(glm.formula),
                          family = 'binomial', maxit = 100) %>% summary()
    Acarrier.model.df <- rbind(Acarrier.model.df, c(Acarrier.model$coefficients[2,1], 
                                                    Acarrier.model$coefficients[,dim(Acarrier.model$coefficients)[2]]))
  }
  colnames(Acarrier.model.df) <- c('allele.COEF.CARRIER', c('Incercept', 'allele', Acarrier.model$coefficients[-c(1,2),] %>% row.names()) %>%
                                     paste('.CARRIER.pval', sep = ''))
  Acarrier.model.df <- data.frame(allele=alleles.freq, Acarrier.model.df)
  
  # Merge 
  glm.data <- merge(Acarrier.model.df, Afreq.model.df, by = 'allele')
  
  # Return
  return(glm.data)
  
}



############### HLA GLM ANALYSIS ################

# Get cases and controls
cases.ids <- covars.df$sample.id.file[which(covars.df$pheno ==1)]
controls.ids <- covars.df$sample.id.file[which(covars.df$pheno ==0)]
data.cases <- HLA.df %>% filter(sample.id.file %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id.file %in% controls.ids)

# Create loci and index
loci <- c("A","B","C","DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
idx <- 1

# Control for allele
controlAllele.df <- controlAllele(settings, loci, HLA.df)

# Iterate over loci for univariate analysis
models.df <- data.frame()

# For each locus 
for (locus in loci){
  
  # Print 
  print(paste0("Current loci: ", locus))
  
  # Subset based on locus
  allele1 <- paste(locus, '.1',sep = '')
  allele2 <- paste(locus,'.2', sep =  '')
  data.locus <- HLA.df[,c("sample.id.file",allele1, allele2)]
  
  # Compute allele frequencies and counts, and carrier frequencies and counts
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  totalCases <-unique(ACFREQ.cases$alleleTotalCase); carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  totalControls <-unique(ACFREQ.controls$alleleTotalControl); carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean 
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df)); ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df));
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
  
  # Parse one hot encoding and merge
  OHE.alleleFreq.data <- alleleFreqOHE(data.locus)
  OHE.alleleFreq.data <- merge(as.data.frame(OHE.alleleFreq.data), 
                               covars.df[c('pheno', 'sample.id.file')], by.x ='sample.id.file', by.y = 'sample.id.file') %>% 
    merge(controlAllele.df, by = 'sample.id.file')
  OHE.carrierFreq.data <- carrierFreqOHE(data.locus)
  OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), 
                                covars.df[c('pheno', 'sample.id.file')], by.x ='sample.id.file', by.y = 'sample.id.file') %>% 
    merge(controlAllele.df, by = 'sample.id.file')
  
  # Run logistic regression model 
  glm.data <- runLogisticRegression(locus, OHE.alleleFreq.data, OHE.carrierFreq.data, covars.df)
  
  # Create dataframes 
  HLA.GLM_alleles.df <-merge(glm.data[,c(1,which(grepl('ALLELE', colnames(glm.data))))],
                             ACFREQ.df[,c(which(grepl(paste(c('A0','A1','A2', 'allele'), collapse = '|'), colnames(ACFREQ.df))))],
                             by = 'allele')
  HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                              ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                              by = 'allele')
  
  # Write to excel output
  if (settings$dataset == "Stanford"){
    write.xlsx(x = HLA.GLM_alleles.df, file = paste(settings$folder$HLA_Output_GLM, 'HLA_AnalysisAlleles','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
    write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$folder$HLA_Output_GLM, 'HLA_AnalysisCarriers','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
  } else if (settings$dataset == "UKB"){
    write.xlsx(x = HLA.GLM_alleles.df, file = paste(settings$folder$HLA_Output_GLM, 'HLA_AnalysisAlleles_UKB','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
    write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$folder$HLA_Output_GLM, 'HLA_AnalysisCarriers_UKB','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
  }
}
