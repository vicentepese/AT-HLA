## ---------------------------
##
## Name: HLA_AA_glm.R
##
## Desciption: Given a dataset of HLA calls, it tests the significance
##              of aminoacid-position pairs based on each allele
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

########## IMPORT #########

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Set values to 0
settings$controlAlleles <- c()
settings$excludePosSubjsByAllele <- c()
settings$allele2Remove <- c()

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls and covariates
HLA.df <- read.csv(settings$file$HLA_df)
covars.df <- read.csv(settings$file$HLA_covars_filt)
covars.df$pheno <- covars.df$pheno -1
covars.df$sex <- covars.df$sex -1

# Add phenotype to HLA datafra,e 
HLA.df <- merge(HLA.df, covars.df[,c("sample.id.file", "pheno")], by = "sample.id.file")

# Import amino acid alignment 
AA_alignment <- read.table(settings$file$AA_alignment, header = TRUE, sep = ',')

# Remove negative positions 
full_prot = FALSE
if (!full_prot){
  AA_alignment[,3] <- AA_alignment %>% apply(MARGIN = 1, function(x) x[3] %>% substr(start = 30, stop = nchar(x[3])))
}

######### CONTROL FOR AMINO ACID ##########

controlAA = function(settings, data, AA_alignment){
  
  # Get locus, position and AA to control 
  AA2control <- settings$AA2control %>% strsplit('_') %>% unlist()
  L2control <- AA2control[1]; posControl <- AA2control[2]; AAcontrol <- AA2control[3]
  c(L2controlA1, L2controlA2) %<-% c(paste0(L2control,'.1'), paste0(L2control,'.2'))
  
  # Get sequences
  AA2control_locus <- AA_alignment %>% filter(locus == AA2control[1])
  
  # Get alleles 
  allelesAA.OG <- AA2control_locus %>% apply(MARGIN = 1, function(x, posControl, AAcontrol)
    if (substr(x[3],posControl,posControl) == AAcontrol) {return(x[2])}, posControl, AAcontrol) %>% unlist()
  allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                         .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist()
  
  # Get data and count cases and controls
  data.AApos <- data %>% filter(get(L2controlA1) %in% allelesAA| get(L2controlA2) %in% allelesAA)
  data.AAneg <- data %>% filter(sample.id.file %notin% data.AApos$sample.if.file)
  
  # Create dataframe with presence of Amino Acid
  data.AAControl <- data.frame(sample.id.file = c(data.AApos$sample.id.file, data.AAneg$sample.id.file),
                               pos = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))))
  colnames(data.AAControl) <- c('sample.id.file', settings$AA2control)
  
  return(data.AAControl)
}

######### DATA TO ONE HOT ENCODING #########

data2OHE = function(settings, data, allelesAA, covars.df, AA_alignment){
  
  # Get data with AAs
  data.AApos <- data %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
  data.AAneg <- data %>% filter(sample.id.file %notin% data.AApos$sample.id.file)
  
  # Create dataframe with OHE
  data.AA <- data.frame(sample.id.file = c(data.AApos$sample.id.file, data.AAneg$sample.id.file), AA = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))),
                        pheno = c(data.AApos$pheno, data.AAneg$pheno))
  
  # Merge with Principal components
  data.AA <- merge(data.AA, covars.df[,-3], by.x = "sample.id.file", by.y = "sample.id.file")
  
  # Merge with presence of controlled amino acid
  if (!settings$AA2control %>% is_empty()){
    data.AAcontrol <- controlAA(settings, data, AA_alignment)
    data.AA <- merge(data.AA, data.AAcontrol, by = "sample.id.file")
  }
  
  # Return 
  return(data.AA)
}

######### RUN GENERALIZED LINEAR MODEL #########

runGLM = function(settings, data.AA){
  
  # GLM formula 
  AA.name <- colnames(data.AA)[2]
  if (length(settings$AA2control ) > 0){
    control.AA <- paste(' ', settings$AA2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
  } else{
    control.AA <- ''
  }
  glm.formula <- paste0("pheno ~ ", AA.name, " + PC1 + PC2 + PC3 ", control.AA)
  
  # Run GLM 
  AA.model <- glm(data = data.AA, formula = as.formula(glm.formula), family = 'binomial', maxit = 100) 
  modelMat <- model.matrix(AA.model)
  
  # Check if the model matrix has full and return coefficients, otherwise return NULL
  rank <- qr(modelMat)$rank
  if (rank == ncol(modelMat)){
    return(AA.model %>% summary())
  } else {
    return(NULL)
  }
}

########## AMINO ACID ANALYSIS ##########

# Get loci, number of cases and controls
loci <- unique(AA_alignment$locus)
c(Ncases, Ncontrols) %<-% c(HLA.df$pheno %>% table %>%.['1'], HLA.df$pheno %>% table %>%.['0'])

# Initialize loop
locus <- c(); AA.eval <- c();  pos.eval <- c(); AA.estimate <- c()
intercept.PVAL <- c(); AA.pval <- c(); PC1.pval <-c(); PC2.pval <- c(); PC3.pval <- c(); AA2control.pval <- c();
alleles <- c(); 

# Initialize loop for map variables
colnameMapPos <- c()
mapPos.df <- matrix(nrow = nrow(HLA.df))

# For each locus 
for (L in loci){
  
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
    AAs <- AAs[which(AAs != '' & AAs!= '*')]
    
    # If more than one, count 
    if (length(AAs) >1){
      
      for (AA in AAs){
        
        # Skip controlled alleles
        if (paste(L,pos,AA, sep='_') == settings$AA2control[1]){
          next
        }
        
        # Get alleles with AAs
        allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
        allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                               .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
        
        # Crete one hot encoding dataframes
        data.AA <- data2OHE(settings, HLA.df, allelesAA, covars.df, AA_alignment); colnames(data.AA)[2] <- paste(L,pos,AA, sep = '_')
        
        # Get map of positions
        mapPos.df <- cbind(mapPos.df, data.AA[,c(paste(L,pos,AA, sep='_'))])
        colnameMapPos <- c(colnameMapPos, paste(L,pos,AA, sep='_'))
        
        # Run generalized linear model 
        model.AA <- runGLM(settings, data.AA)
        
        # If rank is not equal to model matrix (dummy variable trap) skip
        if (is.null(model.AA)){
          next 
        }
        
        # Else, get coefficients
        coefs <- model.AA$coefficients; pvalDim <- dim(coefs)[2]
        
        # Append to vectors
        locus <- c(locus, L); AA.eval <- c(AA.eval, AA); pos.eval <- c(pos.eval, pos);
        AA.estimate <- c(AA.estimate, coefs[2,1]); intercept.PVAL <- c(intercept.PVAL, coefs[1, pvalDim])
        AA.pval <- c(AA.pval, coefs[2, pvalDim]); PC1.pval <- c(PC1.pval, coefs[3, pvalDim]);
        PC2.pval <- c(PC2.pval, coefs[4, pvalDim]); PC3.pval <- c(PC3.pval, coefs[5, pvalDim]); 
        alleles <- c(alleles, paste(allelesAA.OG, collapse = '; '))
        if (dim(coefs)[1] == 6){
          AA2control.pval <- c(AA2control.pval, coefs[6, pvalDim]);
        }
        
      }
    }
  }
}

# Create dataframe
AA.analysis.results <- data.frame(locus= locus, AA = AA.eval, pos = pos.eval, AA.estimate = AA.estimate, intercept.PVAL = intercept.PVAL, 
                                  AA.pval = AA.pval, AA.pvalBONF= AA.pval*length(AA.pval), PC1.pval = PC1.pval, PC2.pval = PC2.pval, PC3.pval = PC3.pval, AA2control.pval = AA2control.pval,
                                  alleles = alleles)
colnames(AA.analysis.results)[ncol(AA.analysis.results)-1] = settings$AA2control

# Write 
write.csv(x = AA.analysis.results, file = paste0(settings$folder$HLA_Output_AA_GLM, "HLA_AA_GLM.csv"))

########## EXTRA FUNCTIONS ############

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




