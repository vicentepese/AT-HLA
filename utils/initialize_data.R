initialize_data = function(settings){

# Import HLA calls, covariates 
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Check pheno
covars.df <- phenoCheck(covars.df)

# If list of matched controls provided, filter
if (!settings$file$matched_controls == ""){
  
  # Verbose
  if (settings$verbose) cat("Parsing matched controls. \n")
  
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
  if (settings$verbose) cat(paste("Ethnicities included:", paste(settings$ethnicity, sep =" "), "\n", sep = " "))

}

# Exclude allele
allele2exclude <- settings$allele2exclude %>% unlist()
if (!allele2exclude %>% is_empty()){
  
  # Verbose 
  if (settings$verbose) cat("Excluded alleles: \n")
  
  for (allele in allele2exclude){

    # Parse locus and allele
    locus <- allele %>% strsplit("\\*") %>% unlist() %>% head(n=1)
    A <- allele %>% strsplit("\\*") %>% unlist() %>% tail(n=1)
    
    # Check allele
    alleleCheck(HLA.df, locus, A)
    
    # Filter HLA calls 
    HLA.df <- HLA.df[which(HLA.df[,paste0(locus,".1")] != A & HLA.df[,paste0(locus,".2")] != A),]
    
    # Verbose
    if (settings$verbose) cat(paste0("\t", allele2exclude,"\n"))
    
  }
}

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

return(list("HLA.df" = HLA.df, "covars.df"=covars.df, "probs.df"=probs.df))
}

