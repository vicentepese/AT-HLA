
# Import libraries
library(jsonlite,warn.conflicts = F)
library(plyr, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse))
library(readr,warn.conflicts = F)
library(data.table, warn.conflicts = F)
library(xlsx,warn.conflicts = F)
library(utils,warn.conflicts = f)
library(ggplot2)

# Set working directory
setwd("~/Documents/HLA_association_pipeline")

# Read settings and import binders 
settings <- jsonlite::read_json("settings.json")
dr7_binders.df <- read_csv("Resources/dr7_binders.csv") 
dr4_binders.df <- read_csv("Resources/dr4_binders.csv") 

# Create command
`%notin%` <- Negate(`%in%`)

# Merge datasets 
binders.df <- dr7_binders.df %>% 
  rbind(dr4_binders.df %>% filter(position %notin% dr7_binders.df$position))

# Create long dataset
binders.long <- pivot_longer(binders.df, cols = 2:18, names_to="DRB1", values_to="rank")
bt <- c()
for (i in 1:nrow(binders.long)){
  if (binders.long[i,3] < 1){
    bt <- c(bt, "Strong binder")
  } else if (binders.long[i,3] < 3){
    bt <- c(bt, "Medium binder")
  } else{
    bt <- c(bt, "Weak binder")
  }
}
binders.long$Binder_type = bt

# Create strong binders dataset and reorder factors
sb <- binders.long %>% filter(bt == "Strong binder" | bt == "Medium binder")
sb$DRB1 <- as.factor(sb$DRB1)
sb$DRB1 <- factor(sb$DRB1, levels = rev(c("07:01", "04:02","04:01", "04:04", "04:05", "04:07", 
                                      "01:01", "03:01", "08:01", "09:01", "11:02", "11:04", "13:01", "13:02","14:01",
                                      "15:01", "16:01")))

binders.long$Binder_type <- factor(binders.long$Binder_type,
                                      levels = c("Strong binder","Medium binder", "Weak binder"))


######################## PLOT #########################

# Bold function
colorado <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

# Main plot
plt <- ggplot(sb %>% filter(Binder_type == "Medium binder")) + 
  geom_point(aes(x=position, y=DRB1, color = Binder_type, cex=1), size = 2.5, shape = 17) + 
  geom_point(data = sb %>% filter(Binder_type == "Strong binder"),
             aes(x=position, y=DRB1, color = Binder_type, cex=1), size = 2.5, shape = 17) +
  xlim(0,520) 
  

# Aesthetics
plt <- plt +
  scale_colour_discrete("Binding affinity") +
  ggtitle("Binding affinities of LGI1 cores to DRB1 alleles associated with anti-LGI1 encephalitis \n and other common alleles") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Binding core position") +
  ylab("DRB1") + 
  theme(axis.text.y=element_text(face=colorado(sb$DRB1, c("07:01", "04:02")))) + 
  coord_fixed(ratio=12)
plt
ggsave("binding_affinities_xlim.png",dpi=600)


# Add LGI1 sections
plt + geom_vline(xintercept = 34) + geom_vline(xintercept = 72) +
  geom_vline(xintercept = 92) + geom_vline(xintercept = 113)  + # LRR1
  geom_vline(xintercept = 116) + geom_vline(xintercept = 137) + # LRR2
  geom_vline(xintercept = 140) + geom_vline(xintercept = 161) + # LRR3
  geom_vline(xintercept = 173) + geom_vline(xintercept = 223) + # LRRCT
  geom_vline(xintercept = 225) + geom_vline(xintercept = 267) + # EAR 1
  geom_vline(xintercept = 271) + geom_vline(xintercept = 313) + # Ear 2
  geom_vline(xintercept = 317) + geom_vline(xintercept = 364) + # EAR 3
  geom_vline(xintercept = 366) + geom_vline(xintercept = 415) + # EAR 4
  geom_vline(xintercept = 419) + geom_vline(xintercept = 462) + # EAR5
  geom_vline(xintercept = 464) + geom_vline(xintercept = 506) + # EAR 6
  geom_vline(xintercept = 512) + geom_vline(xintercept = 552)
ggsave("binding_w_pos.png",dpi=600)
``