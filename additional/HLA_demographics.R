# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(xlsx)
library(plyr)
library(epitools)

########### INITIALIZATION ########### 

# Set working directory
setwd("~/Documents/HLA_association_pipeline")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls, covariates 
HLA.df <- read.csv(settings$file$HLA_Data)
covars.df <- read.csv(settings$file$covars)
probs.df <- read.csv(settings$file$probs)

# Read options
prob_thr <- settings$prob_thr
freq_thr <- settings$freq_thr*100
alleles2control <- settings$allele2exclude %>% unlist()

# Filter out the alleles to exclude 
for (allele in alleles2control){
  
  # Parse locus and allele
  locus <- allele %>% strsplit("\\*") %>% unlist() %>% head(n=1)
  A <- allele %>% strsplit("\\*") %>% unlist() %>% tail(n=1)
  
  # Filter HLA calls 
  HLA.df <- HLA.df[which(HLA.df[,paste0(locus,".1")] != A & HLA.df[,paste0(locus,".2")] != A),]
  
}

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

# Import demographics
demo.df <- read.table("Resources/demographics.csv", header = TRUE, sep = ",")

############## SEX ANALYSIS ###############

# Get HLA calls with sex
HLA.df.sex <- HLA.df %>% merge(demo.df, by = "sample.id")
HLA.df.sex <- HLA.df.sex %>% merge(covars.df[,c("sample.id", "pheno")], by = "sample.id")

# Get hetero, homo and negative dataset of allele to study 
a2study <- settings$allele2control %>% unlist() 
locus <- a2study %>% strsplit("\\*") %>% unlist() %>%  head(n=1); locus.ids <- paste0(rep(locus,2), c(".1",".2"))
allele <- a2study %>% strsplit("\\*") %>% unlist() %>%  tail(n=1)
HLA.het.sex <- HLA.df.sex[xor(x = HLA.df.sex[locus.ids[1]] == allele, y = HLA.df.sex[locus.ids[2]] == allele),]
HLA.homo.sex <-  HLA.df.sex %>% filter(get(locus.ids[1]) == allele & get(locus.ids[2]) == allele)
HLA.neg.sex <- HLA.df.sex %>% filter(sample.id %notin% HLA.homo.sex$sample.id & sample.id %notin% HLA.het.sex$sample.id)

# Compute contingency table
het.male <- HLA.het.sex$Sex %>% table() %>%.["M"]; het.fem <- HLA.het.sex$Sex %>% table() %>%.["F"]
homo.male <- HLA.homo.sex$Sex %>% table() %>% .["M"]; homo.fem <- HLA.homo.sex$Sex %>% table() %>% .["F"]
neg.male <- HLA.neg.sex$Sex %>% table() %>%.["M"]; neg.fem <- HLA.neg.sex$Sex %>% table() %>%.["F"]
totalMales <- het.male + homo.male +  neg.male; totalFemales <- het.fem + homo.fem + neg.fem
cont.table.homo <- matrix(c(het.male, het.fem, homo.male, homo.fem), nrow = 2)
cont.table.neg <- matrix(c(het.male, het.fem, neg.male, neg.fem), nrow = 2)

# Compute Chi2
chisq.res.homo <- chisq.test(cont.table.homo); OR.res.homo <- oddsratio.fisher(cont.table.homo)
chisq.res.neg <- chisq.test(cont.table.neg); OR.res.neg <- oddsratio.fisher(cont.table.neg)

# Create dataframe 
sex.df <- data.frame(allele1 = c(a2study, a2study, "X"), allele2= c(a2study, "X", "X"), 
                     Males = c(homo.male, het.male, neg.male),
                     FreqMales = c(homo.male, het.male, neg.male)/totalMales*100, 
                     Females = c(homo.fem, het.fem, neg.fem), 
                     FreqFemales = c(homo.fem, het.fem, neg.fem)/totalFemales*100,
                     Chi2.PVAL.HET = c(NA, chisq.res.neg$p.value, "Ref"),
                     OR = c(NA, paste0(OR.res.neg$measure[2,1]), "Ref"),
                     OR_LI = c(NA, paste0(OR.res.neg$measure[2,2]), "Ref"),
                     OR_UI = c(NA, paste0(OR.res.neg$measure[2,3]), "Ref"),
                     Chi2.PVAL.HOMO = c(chisq.res.homo$p.value, "Ref", NA), 
                     OR = c(OR.res.homo$measure[2,1], "Ref", NA),
                     OR_LI = c(OR.res.homo$measure[2,2], "Ref", NA),
                     OR_UI = c(OR.res.homo$measure[2,3], "Ref", NA))
# Write output
write.csv(sex.df, paste0(settings$Output$Utils, "Zygosity_demographics.csv"), row.names = FALSE)

############## SEX ANALYSIS in heterozygous ###############

# Parse other alleles in the heterozygous group, and remove allele 2 study
allelesX <- HLA.het.sex[,c(locus.ids)] %>% unlist() %>% unique() 
allelesX <- allelesX[-which(allelesX==allele)]

# For the rest of alleles, compute fishers
hetero.df <- data.frame()
for (allele in allelesX){
  
  # Get count 
  HLA.het.allele <- HLA.het.sex %>% filter(get(locus.ids[1]) == allele | get(locus.ids[2]) == allele)
  
  # Count males and females
  aMale <- HLA.het.allele$Sex %>% table() %>% .["M"]; aFemale <-HLA.het.allele$Sex %>% table() %>% .["F"]
  
  # Compute contingency table 
  cont.table <- matrix(c(aMale, aFemale, totalMales-aMale, totalFemales-aFemale),2)
  cont.table[is.na(cont.table)] <- 0
  
  # Fishers 
  fisher.res <- fisher.test(cont.table)
  
  # Append 
  hetero.df <- hetero.df %>% rbind(c(allele,aMale, aFemale, fisher.res$p.value, 
                                   fisher.res$estimate, fisher.res$conf.int[1], fisher.res$conf.int[2]))
  
}

# Change colnames
colnames(hetero.df) <- c("allele", "NMale", "NFemale","Fisher.PVAL","OR", "LCI","UCI")

# Write output
write.csv(hetero.df, paste0(settings$Output$Utils, "Heterozygosity_demographics.csv"), row.names = FALSE)


