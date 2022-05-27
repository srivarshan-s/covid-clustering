#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")



#################### IMPORT PACKAGES ###########################
library("dplyr")
library("tidyverse")
library("readr")



#################### MAIN CODE #################################

# Read data
df <- read_tsv("data/ECG200_TRAIN.tsv", col_names = FALSE)
print(df)
