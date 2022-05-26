# Environment Setup
setwd("~/Documents/covid-clustering")

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
        "stringency_index", "new_tests_per_thousand")
# print(df)

# Filter the countries
countries <- read_csv("data/countries.csv")
df2 <- merge(df, countries, by.x = "iso_code", by.y = "country.code")
print(df2)