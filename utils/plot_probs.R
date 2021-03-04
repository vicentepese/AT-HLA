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
library(gridExtra)

########### INITIALIZATION ############

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


############## PLOT PROBABILITIES #############

## Plot probabilities distributions
pl <- vector('list', ncol(probs.df)-3)
idx <- 1
for (i in 4:ncol(probs.df)){
  
  pl[[i-3]] <- local({
    i <- i
    p1 <- ggplot(probs.df, aes(get(colnames(probs.df)[i]))) +
      geom_histogram() +
      xlab(colnames(probs.df)[i]) + xlim(0,1)
    print(p1)
  })
}
do.call(grid.arrange, pl)

##  Plot PCA probabilities

# Merge PCA with probs 
data_plot <- merge(covars.df, probs.df, by = c("FID", "IID","sample.id"))

# Crate array of probs.IDs
loci <- c("A","B","C","DPA1","DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
prob.ids <- paste0(rep("prob.", length(loci)), loci)

# Plot PCAs with probs
pl <- vector('list', length(loci))
idx <- 1
for (i in prob.ids){
  
  pl[[idx]] <- local({
    i <- i
    p1 <- ggplot(data_plot, aes(PC1, PC2, color = get(i))) + geom_point(alpha = 0.6) +
      scale_colour_gradientn(limits = c(0,1), colors =c("navyblue", "darkmagenta", "darkorange1"), oob = scales::squish) + 
      labs(color = i)
    print(p1)
  })
  idx <- idx+1
}
do.call(grid.arrange, pl)
g <- arrangeGrob(pl) #generates g
ggsave(file="pca_probs.pdf",g)