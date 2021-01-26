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
library(TreeBH)

########## IMPORT #########

# Set working directory
setwd("~/Documents/HLA_association_pipeline")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Import data 
HLA.df  <- read.csv(settings$file$HLA_Data)

# Impport covariates
covars.df <- read.csv(settings$file$covars)

# Not-in command
`%notin%` <- Negate(`%in%`)

########### CONTROL ALLELES ##########

controlAllele = function(as2control, HLA.df){
  
  # Get unique loci 
  lociControl = lapply(as2control, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist() %>% unique() 
  
  # For each allele 
  alleleControl.df = data.frame(sample.id.file = HLA.df$sample.id)
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
  return(alleleControl.df)
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


######### GLM ###########

runLogisticRegression = function(locus, OHE.carrierFreq.data, covars.df, as2control = NULL){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  covars.df$sample.id <- as.factor(covars.df$sample.id)
  alleles.freq <- colnames(OHE.carrierFreq.data)[-c(1,(ncol(OHE.carrierFreq.data)-length(as2control)):ncol(OHE.carrierFreq.data))]
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, covars.df[,-3], by = 'sample.id')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = locus)]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  for (allele in alleles.freq){
    if (!is.null(as2control)){
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
  glm.data <- Acarrier.model.df
  
  # Return
  return(glm.data)
  
}

fitGLM = function(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control = NULL){
  
  # Compute allele frequencies and counts, and carrier frequencies and counts
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
  
  # Control for allele
  controlAllele.df <- controlAllele(as2control, HLA.df)
  
  # Subset locus
  allele1 <- paste(locus, '.1',sep = '')
  allele2 <- paste(locus,'.2', sep =  '')
  data.locus <- HLA.df[,c('sample.id', allele1, allele2)]

  # Parse one hot encoding and merge
  OHE.carrierFreq.data <- carrierFreqOHE(data.locus); if ('V1' %in% colnames(OHE.carrierFreq.data)){OHE.carrierFreq.data$V1 <- NULL}
  OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), 
                                covars.df[c('pheno', 'sample.id')], by.x ='sample.id', by.y = 'sample.id') %>% 
    merge(controlAllele.df, by = 'sample.id')
  
  # Remove subjects thar are not controls or cases
  OHE.carrierFreq.data <- OHE.carrierFreq.data %>% filter(pheno != -9)
  
  # Run logistic regression model 
  glm.data <- runLogisticRegression(locus, OHE.carrierFreq.data, covars.df, as2control)
  
  # Create dataframes 
  HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                              ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                              by = 'allele')
  
  return(HLA.GLM_carriers.df)
}


########## HLA ANALYSIS ###########

# Get cases and controls
cases.ids <- covars.df$sample.id[which(covars.df$pheno ==1)]
controls.ids <- covars.df$sample.id[which(covars.df$pheno ==0)]
data.cases <- HLA.df %>% filter(sample.id %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id %in% controls.ids)

# Initialize while lopp
pval <- 0; as2control <- c(); signAlleles <- list(); 

# HLA Loci
loci <- c('A','B','C','DQA1', 'DQB1', 'DPB1', 'DRB1','DRB3','DRB4','DRB5')

# While signifiant alleles
idx <- 1; iter <- 1
while (pval < 0.05){
  
  # Initialize loop over locus
  HLA.CHI_carriers.list <- list(); HLA.GLM_carriers.list <- list(); 
  pvalTotal <- c()
  
  # Compute Chi2 and fit glm for each locus
  for (locus in loci){
    
    # Fit GLM
    HLA.GLM_carriers.df <- fitGLM(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control)

    # Add to final list
    HLA.GLM_carriers.list[[locus]] <- HLA.GLM_carriers.df 

    # Add p-values
    for (A in HLA.GLM_carriers.df$allele){
      pvalA <- HLA.GLM_carriers.df$allele.CARRIER.pval[which(HLA.GLM_carriers.df$allele == A)]
      pvalTotal[paste(locus, '*', A, sep = '')] <- pvalA
    }
    
    # Get minimum allele value 
    pval <- pvalTotal[which(pvalTotal == min(pvalTotal))][1]
    pvalMin <- pval %>% names(); pvalMinLocus <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[1]
    pvalMinAllele <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Add to outputs
    as2control <- c(as2control, pvalMin)
    preOut <- HLA.GLM_carriers.list[[ pvalMinLocus]] %>% filter(allele == pvalMinAllele)
    preOut$allele <- pvalMin  ; signAlleles[[idx]] <- preOut
    
    # Update
    iter <- iter + 1
    idx <- idx +1
  }
}

# Format output
allele <- as.data.frame(signAlleles[[1]]); 
for (idx in 2:length(signAlleles)){
  allele <- rbind.fill(allele, signAlleles[[idx]] %>% as.data.frame())
}

# Write
for (idx in 1:length(HLA.GLM_carriers.list)){
  HLA.GLM_carriers.df <- HLA.GLM_carriers.list[[idx]]; locus <- HLA.GLM_carriers.list %>% names() %>% .[idx]
  write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$folder$HLA_Output, 'HLA_GLM_Carriers','.xlsx', sep = ''), sheetName = locus, 
             col.names = TRUE, row.names = FALSE, append = TRUE)
}
write.xlsx(x = allele, file = paste(settings$folder$HLA_Output, 'HLA_GLM_Carriers','.xlsx', sep = ''), sheetName = 'Significant_alleles', 
           col.names = TRUE, row.names = FALSE, append = TRUE)
