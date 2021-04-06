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

#Read demographics, merge with covars, and filter HLA
demo.df <- read.table("Resources/demographics.csv", sep = ",", header = TRUE)
demo.df$Age[is.na(demo.df$Age)] <- 0; demo.df$OnsetAge[is.na(demo.df$OnsetAge)] <- 0
demo.df$Age <- demo.df$Age+demo.df$OnsetAge; demo.df$Age[demo.df$Age==0] <- NA
covars.df <- covars.df %>% merge(demo.df[,c("sample.id", "Age", "OnsetAge", "Sex")], by = "sample.id")

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

############### T-TEST ############

# Parse allele of interest
a2control <- settings$allele2control %>% unlist() %>% head(n=1)
locus <- a2control %>% strsplit("\\*") %>% unlist() %>% head(n=1); locus.ids <- paste0(rep(locus,2), c(".1",".2"))
allele <- a2control %>% strsplit("\\*") %>% unlist() %>% tail(n=1)

# Parse heterozygous and homozygous HLA subjects
HLA.hetero <- HLA.df[xor(HLA.df[,locus.ids[1]] == allele, HLA.df[,locus.ids[2]] == allele),]
HLA.homo <- HLA.df %>% filter(get(locus.ids[1]) == allele & get(locus.ids[2]) == allele)
HLA.neg <- HLA.df %>% filter(get(locus.ids[1]) != allele & get(locus.ids[2]) != allele)

# Merge and add age 
HLA.hetero <- HLA.hetero %>% merge(demo.df[,c("sample.id", "Sex")])
HLA.homo <- HLA.homo %>% merge(demo.df[,c("sample.id", "Sex")])
HLA.neg <- HLA.neg %>% merge(demo.df[,c("sample.id", "Sex")])

# Parse by groups
HLA.G1 <- HLA.homo; G1.male <- HLA.G1$Sex %>% table() %>% .["M"]; G1.female <- HLA.G1$Sex %>% table() %>% .["F"];
HLA.G2 <- HLA.hetero %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02");  G2.male <- HLA.G2$Sex %>% table() %>% .["M"]; G2.female <- HLA.G2$Sex %>% table() %>% .["F"];
HLA.G3 <- HLA.neg %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02"); G3.male <- HLA.G3$Sex %>% table() %>% .["M"]; G3.female <- HLA.G3$Sex %>% table() %>% .["F"];
HLA.G4 <- HLA.hetero %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02"); G4.male <- HLA.G4$Sex %>% table() %>% .["M"]; G4.female <- HLA.G4$Sex %>% table() %>% .["F"];
HLA.G5<- HLA.neg %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02"); G5.male <- HLA.G5$Sex %>% table() %>% .["M"]; G5.female <- HLA.G5$Sex %>% table() %>% .["F"];

# # Compute t-test 
# hetero.test <- t.test(HLA.hetero$Sex[!is.na(HLA.hetero$Sex)], HLA.neg$Sex[!is.na(HLA.neg$Sex)])
# homo.test <- t.test(HLA.hetero$Sex[!is.na(HLA.hetero$Sex)], HLA.homo$Sex[!is.na(HLA.homo$Sex)])

# Compute t-test
G3.test <- fisher.test(matrix(c(G3.male, G4.male, G3.female, G4.female), nrow = 2))
G2.test <- fisher.test(matrix(c(G2.male, G4.male, G2.female, G4.female), nrow = 2))
G5.test <- fisher.test(matrix(c(G5.male, G4.male, G5.female, G4.female), nrow = 2))
G1.test <- fisher.test(matrix(c(G1.male, G4.male, G1.female, G4.female), nrow = 2))

# Create dataframe 
Sex.df <- data.frame(allele1 = c("DR7-", "DR7+", "DR7-","DR7+","DR7+"), allele2= c("DR4+", "DR4+","DR4-","DR7+","DR4-"),
                     NMale = c(G3.male, G2.male, G5.male, G1.male, G4.male), 
                     FreqMale = c(G3.male, G2.male, G5.male, G1.male, G4.male)/nrow(HLA.df),
                     NFemale = c(G3.female, G2.female, G5.female, G1.female, G4.female),
                     FreqFemale = c(G3.female, G2.female, G5.female, G1.female, G4.female)/nrow(HLA.df),
                     fisher.pval = c(G3.test$p.value, G2.test$p.value, G5.test$p.value, G1.test$p.value, "Ref"),
                     OR= c(G3.test$estimate, G2.test$estimate, G5.test$estimate, G1.test$estimate, "Ref"), 
                     LCI = c(G3.test$conf.int[1], G2.test$conf.int[1], G5.test$conf.int[1], G1.test$conf.int[1], "Ref"),
                     UCI = c(G3.test$conf.int[2], G2.test$conf.int[2], G5.test$conf.int[2], G1.test$conf.int[2], "Ref")
)


# Write 
write.xlsx(x = Sex.df, file = "Outputs/Utils/Age_T_Sex.xlsx", sheet = "Sheet1")

