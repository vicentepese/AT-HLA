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
demo.df <- read.table("Resources/mrs_scores.csv", sep = ",", header = TRUE)
demo.df$FBDS <- as.numeric(demo.df$FBDS); demo.df$mrs.onset <- as.numeric(demo.df$mrs.onset)
demo.df$mrs.max <- as.numeric(demo.df$mrs.max); demo.df$mrs.last <- as.numeric(demo.df$mrs.last)
covars.df <- covars.df %>% merge(demo.df %>% select(-one_of("FID", "IID")), by = "sample.id")

# Parse HLA calls for which there is a phenotype 
HLA.df <- HLA.df %>% filter(sample.id %in% covars.df$sample.id)

############### T-TEST ############

# Parse allele of interest
a2control <- settings$allele2control %>% unlist() %>% head(n=1)
locus <- a2control %>% strsplit("\\*") %>% unlist() %>% head(n=1); locus.ids <- paste0(rep(locus,2), c(".1",".2"))
allele <- a2control %>% strsplit("\\*") %>% unlist() %>% tail(n=1)

# Colnames
cols = colnames(demo.df)[4:ncol(demo.df)]

# Parse heterozygous and homozygous HLA subjects
HLA.hetero <- HLA.df[xor(HLA.df[,locus.ids[1]] == allele, HLA.df[,locus.ids[2]] == allele),]
HLA.homo <- HLA.df %>% filter(get(locus.ids[1]) == allele & get(locus.ids[2]) == allele)
HLA.neg <- HLA.df %>% filter(get(locus.ids[1]) != allele & get(locus.ids[2]) != allele)

# Merge and add age 
HLA.hetero <- HLA.hetero %>% merge(demo.df[,c("sample.id", cols)])
HLA.homo <- HLA.homo %>% merge(demo.df[,c("sample.id", cols)])
HLA.neg <- HLA.neg %>% merge(demo.df[,c("sample.id", cols)])

# Parse by groups
HLA.G1 <- HLA.homo
HLA.G2 <- HLA.hetero %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02")
HLA.G3 <- HLA.neg %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02")
HLA.G4 <- HLA.hetero %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02")
HLA.G5<- HLA.neg %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02")

# Initialize output dataframe
mrs.df <- data.frame(allele1 = c("DR7-", "DR7+", "DR7-","DR7+","DR7+"), allele2= c("DR4+", "DR4+","DR4-","DR7+","DR4-"),
                     Ncases = c(nrow(HLA.G3), nrow(HLA.G2), nrow(HLA.G5), nrow(HLA.G1), nrow(HLA.G4)))
mrsHH.df <- data.frame(allele1 = c(allele, allele, "X"), allele2= c(allele, "X", "X"))

# Compute t-test
for (col in cols){
  
  # Compute means and stds
  means <- c(mean(HLA.G3[,col], na.rm=T), mean(HLA.G2[,col], na.rm=T), mean(HLA.G5[,col], na.rm=T),
            mean(HLA.G1[,col], na.rm=T), mean(HLA.G4[,col], na.rm=T))
  stds <- c(sd(HLA.G3[,col], na.rm=T), sd(HLA.G2[,col], na.rm=T), sd(HLA.G5[,col], na.rm=T),
            sd(HLA.G1[,col], na.rm=T), sd(HLA.G4[,col], na.rm=T))
  
  # Compute t test
  G3.t <- t.test(HLA.G3[,col], HLA.G4[,col]) 
  G2.t <- t.test(HLA.G2[,col], HLA.G4[,col])
  G5.t <- t.test(HLA.G5[,col], HLA.G4[,col]) 
  G1.t <- t.test(HLA.G1[,col], HLA.G4[,col]) 
  
  # Compute vectors 
  t.stat <-  c(G3.t$statistic, G2.t$statistic, G5.t$statistic, G1.t$statistic, "Ref")
  p.vals <-c (G3.t$p.value, G2.t$p.value, G5.t$p.value, G1.t$p.value, "Ref")
  
  # Add columns
  mrs.df <- mrs.df %>% add_column(means); colnames(mrs.df)[ncol(mrs.df)] <- paste0(col,"_mean")
  mrs.df <- mrs.df %>% add_column(stds); colnames(mrs.df)[ncol(mrs.df)] <- paste0(col,"_std")
  mrs.df <- mrs.df %>% add_column(t.stat); colnames(mrs.df)[ncol(mrs.df)] <- paste0(col,"_t.stat")
  mrs.df <- mrs.df %>% add_column(p.vals); colnames(mrs.df)[ncol(mrs.df)] <- paste0(col,"_p.val")
  
  # Compute means and stds HH
  means <- c(mean(HLA.homo[,col], na.rm=T), mean(HLA.hetero[,col], na.rm=T), mean(HLA.neg[,col], na.rm=T))
  stds <- c(sd(HLA.homo[,col], na.rm=T), sd(HLA.hetero[,col], na.rm=T), sd(HLA.neg[,col], na.rm=T))
  
  # Compute t test HH
  homo.t <- t.test(HLA.homo[,col], HLA.hetero[,col])
  hetero.t <- t.test(HLA.hetero[,col], HLA.neg[,col])
  
  # Compute vectors HH
  T.value.HOMO = c(homo.t$statistic, "Ref", NA)
  T.PVAL.HOMO = c(homo.t$p.value, "Ref", NA)
  T.value.HET = c(NA, hetero.t$statistic, "Ref")
  T.PVAL.HET = c(NA, hetero.t$p.value, "Ref")
  
  # Add columns
  mrsHH.df <-  mrsHH.df %>% add_column(means); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_mean")
  mrsHH.df <-  mrsHH.df %>% add_column(stds); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_stds")
  mrsHH.df <- mrsHH.df %>% add_column(T.value.HOMO); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_t.stat_HOMO")
  mrsHH.df <- mrsHH.df %>% add_column(T.PVAL.HOMO); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_p.val_HOMO")
  mrsHH.df <- mrsHH.df %>% add_column(T.value.HET); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_t.stat_HET")
  mrsHH.df <- mrsHH.df %>% add_column(T.PVAL.HET); colnames(mrsHH.df)[ncol(mrsHH.df)] <- paste0(col,"_p.val_HET")
  
}

# Write 
write.xlsx(x = mrs.df, file = "Outputs/Utils/mRS_t_test.xlsx", sheet = "Sheet1")
write.xlsx(x = mrsHH.df, file = "Outputs/Utils/mRSHH_t_test.xlsx", sheet = "Sheet1")

