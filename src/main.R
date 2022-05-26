# Environment setup
setwd("~/Documents/covid-clustering")

# Source other scripts
source("src/negatives.R")

# Import packages
library(readr)
library(dplyr)
library(tidyverse)

# Importing the data
df <- read_csv("data/owid_covid_data.csv")

# Print the data
# print(df)

# Select only the required columns
df <- df %>%
    select(
        "iso_code", "continent", "location", "date",
        "new_cases_per_million", "new_deaths_per_million",
        "total_cases_per_million", "total_deaths_per_million",
        "stringency_index", "new_tests_per_thousand"
    )
# print(df)

# Filter the countries
countries <- read_csv("data/countries.csv")
df <- merge(df, countries, by.x = "iso_code", by.y = "country.code")
# print(df)

# Set missing values (na) to 0
df[is.na(df)] <- 0
# print(df)

# Fix negatives in the columns
df <- fix_negatives(fix_data = df, column_to_fix = "new_cases_per_million")
print(df)
