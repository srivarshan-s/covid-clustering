# Environment setup
setwd("~/Documents/chemicals/visualization-tools/")
source("../open_chemicals.R")

# Load packages
packages <- c("palmerpenguins", "ochRe", "GGally", 
              "tourr", "liminal", "tidyverse")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Load data
chemical_labels <- read_csv("../chemical_index.csv")
