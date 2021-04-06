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
HLA.hetero <- HLA.hetero %>% merge(demo.df[,c("sample.id", "Age")])
HLA.homo <- HLA.homo %>% merge(demo.df[,c("sample.id", "Age")])
HLA.neg <- HLA.neg %>% merge(demo.df[,c("sample.id", "Age")])

# Parse by groups
HLA.G1 <- HLA.homo
HLA.G2 <- HLA.hetero %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02")
HLA.G3 <- HLA.neg %>% filter(DRB1.1 == "04:02" | DRB1.2 == "04:02")
HLA.G4 <- HLA.hetero %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02")
HLA.G5<- HLA.neg %>% filter(DRB1.1 != "04:02" & DRB1.2 != "04:02")

# # Compute t-test 
# hetero.test <- t.test(HLA.hetero$Age[!is.na(HLA.hetero$Age)], HLA.neg$Age[!is.na(HLA.neg$Age)])
# homo.test <- t.test(HLA.hetero$Age[!is.na(HLA.hetero$Age)], HLA.homo$Age[!is.na(HLA.homo$Age)])

# Compute t-test
G3.t <- t.test(HLA.G3$Age[!is.na(HLA.G3$Age)], HLA.G4$Age[!is.na(HLA.G4$Age)])
G2.t <- t.test(HLA.G2$Age[!is.na(HLA.G2$Age)], HLA.G4$Age[!is.na(HLA.G4$Age)])
G5.t <- t.test(HLA.G5$Age[!is.na(HLA.G5$Age)], HLA.G4$Age[!is.na(HLA.G4$Age)])
G1.t <- t.test(HLA.G1$Age[!is.na(HLA.G1$Age)], HLA.G4$Age[!is.na(HLA.G4$Age)])

# Create dataframe 
age.df <- data.frame(allele1 = c("DR7-", "DR7+", "DR7-","DR7+","DR7+"), allele2= c("DR4+", "DR4+","DR4-","DR7+","DR4-"),
                     Ncases = c(nrow(HLA.G3), nrow(HLA.G2), nrow(HLA.G5), nrow(HLA.G1), nrow(HLA.G4)), 
                     Freq = c(nrow(HLA.G3), nrow(HLA.G2), nrow(HLA.G5), nrow(HLA.G1), nrow(HLA.G4))/nrow(HLA.df),
                     mean.Age = c(mean(HLA.G3$Age, na.rm = T), mean(HLA.G2$Age, na.rm = T), mean(HLA.G5$Age, na.rm = T), 
                                  mean(HLA.G1$Age, na.rm = T),mean(HLA.G4$Age, na.rm = T)),
                     std.Age = c(sd(HLA.G3$Age, na.rm = T), sd(HLA.G2$Age, na.rm = T), sd(HLA.G5$Age, na.rm = T), 
                                 sd(HLA.G1$Age, na.rm = T),sd(HLA.G4$Age, na.rm = T)),
                     t.stat = c(G3.t$statistic, G2.t$statistic, G5.t$statistic, G1.t$statistic, "Ref"),
                     p.val= c(G3.t$p.value, G2.t$p.value, G5.t$p.value, G1.t$p.value, "Ref")
                     )

age.df <- data.frame(allele1 = c(allele, allele, "X"), allele2= c(allele, "X", "X"),
                     mean.Age = c(mean(HLA.homo$Age, na.rm = T), mean(HLA.hetero$Age, na.rm = T), mean(HLA.neg$Age, na.rm = T)),
                     std.Age = c(sd(HLA.homo$Age, na.rm = T), sd(HLA.hetero$Age, na.rm = T), sd(HLA.neg$Age, na.rm = T)),
                     T.value.HOMO = c(homo.test$statistic, "Ref", NA),
                     T.LI.HOMO = c(homo.test$conf.int[1], "Ref", NA), T.UI.HOMO = c(homo.test$conf.int[2], "Ref", NA),
                     T.PVAL.HOMO = c(homo.test$p.value, "Ref", NA),
                     T.value.HET = c(NA, hetero.test$statistic, "Ref"),
                     T.LI.HET = c(NA, hetero.test$conf.int[1], "Ref"), T.UI.HET = c(NA, hetero.test$conf.int[2], "Ref"),
                     T.PVAL.HOMO = c(NA, hetero.test$p.value, "Ref"))

# Write 
write.xlsx(x = age.df, file = "Outputs/Utils/Age_T_Test.xlsx", sheet = "Sheet1")

