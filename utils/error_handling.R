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
  if (settings$prob_thr < 0 | settings$prob_thr >= 1){
    stop("The probability threshold must be between 0 and 1.")
  }
  if (settings$freq_thr < 0 | settings$prob_thr >=1){
    stop("The frequency threshold must be between 0 and 1.")
  }
  
  # Throw warning 
  if (settings$prob_thr == 0){
    warning("Probability threshold is set to 0. All cases will be included.")
  } 
  if (settings$freq_thr == 0){
    warning("Frequency threshold is set to 0. All alleles will be included.")
  }
  
}

phenoCheck = function(covars.df){

  # Get phenotype
  pheno_labels <- covars.df$pheno %>% unique()

  # Check that phenotype is between 0 and 1, or 1 and 2 
  if (all(pheno_labels %in% c(1,2))){  
    covars.df$pheno <- covars.df$pheno - 1
  } else if (any(pheno_labels %notin% c(0,1))){ 
    stop("Phenotype labels are not 1 and 2, or 0 and 1.")
  }
  
  # Get number of cases and controls
  pheno_count <- covars.df$pheno  %>% table()
  
  # Check that there are cases and controls
  if (pheno_count['0'] == 0){
    stop("No controls in the covariates file.")
  } else if (pheno_count['1'] == 0){
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
  
  # Check locus
  if (!any(grepl(locus, HLA.df %>% colnames()))){
    stop("Locus in allele to exclude/controls not present in the data.")
  }

  # Get alleles
  alleles <- unique(c(HLA.df[,c(paste0(locus,".1"))], HLA.df[,c(paste0(locus,".2"))]))
  
  # Check that allele existst
  if (A %notin% alleles){
    stop("Allele to exclude/control not present in the data.")
  }

}

ethnicityCheck = function(settings, ethnicity.df){
  
  # Check that ethnicities are in the data provided 
  if (any(settings$ethnicity %notin% ethnicity.df$Population)){
    stop("One of the provided ethnicities is not present in the data.")
  }
  
}