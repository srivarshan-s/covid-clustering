#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")



#################### IMPORT PACKAGES ###########################
library("dplyr")
library("tidyverse")
library("readr")



#################### MAIN CODE #################################

# Read train data
df_train <- read_tsv("data/ECG200_TRAIN.tsv", col_names = FALSE)
# print(df_train)

# Read test data
df_test <- read_tsv("data/ECG200_TEST.tsv", col_names = FALSE)
# print(df_test)

# Merging the two dataframes
df <- rbind(df_train, df_test)
# print(df)

# Extracting the two classes
df_class_1 <- filter(df, X1 == "-1")
df_class_2 <- filter(df, X1 == "1")
print(df_class_1)
print(df_class_2)