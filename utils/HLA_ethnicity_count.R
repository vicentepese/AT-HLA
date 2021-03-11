## ---------------------------
##
## Name: HLA_ethnicity_count.R
##
## Desciption: Given a dataset of HLA calls, an allele and
##    an ethnicity, it counts the number of cases/controls per ethnicity.
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
library(xlsx)
library(plyr)