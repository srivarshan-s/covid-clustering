# Environment Setup
setwd("~/Documents/chemicals/")
source("open_chemicals.R")
library(tourr)



# Helper Functions
animate_single_coeff <- function (coeff, guide) {
  chem_df <- data.frame(matrix(nrow = 63, ncol = 12))
  chem_names <- c()
  tour_dependence <- c()
  for (i in 1:11) {
    chem_df[, i] <- as.numeric(chem$fd[[i]]$coefs[coeff,])
    chem_names <- append(chem_names, paste0("C", i, ".", coeff))
    tour_dependence <- append(tour_dependence, c(i))
  }
  names(chem_df) <- append(chem_names, "category")
  chem_df$category <- ifelse(chem$groupd == 1 | chem$groupd == 10,
                             as.character(chem$groupd), "other")
  dev.new()
  if (guide) {
    animate_xy(chem_df[, 1:11], col = chem_df$category, fps = 30,
               tour_path = guided_tour(holes()))
  } else {
    animate_xy(chem_df[, 1:11], col = chem_df$category, fps = 30)
  }
}

animate_multiple_coeff <- function (coeffs, guide) {
  len <- length(coeffs)
  ncol = 1 + 11 * len
  chem_df <- data.frame(matrix(nrow = 63, ncol = ncol))
  chem_names <- c()
  tour_dependence <- c()
  for (j in 1:len) {
    for (i in 1:11) {
      i_1 <- i + 11 * (j - 1)
      chem_df[, i_1] <- as.numeric(chem$fd[[i]]$coefs[coeffs[j],])
      chem_names <- append(chem_names, paste0("C", i, ".", coeffs[j]))
      tour_dependence <- append(tour_dependence, c(i))
    }
  }
  names(chem_df) <- append(chem_names, "category")
  chem_df$category <- ifelse(chem$groupd == 1 | chem$groupd == 10,
                             as.character(chem$groupd), "other")
  dev.new()
  if (guide) {
    animate_xy(chem_df[, 1:ncol-1], col = chem_df$category, fps = 30,
               tour_path = guided_tour(holes()))
  } else {
    animate_xy(chem_df[, 1:ncol-1], col = chem_df$category, fps = 30)
  }
}

animate_all_coeff <- function(num_splines, guide) {
  chem_df <- data.frame(matrix(nrow = 63*num_splines, ncol = 12))
  chem_names <- c()
  tour_dependence <- c()
  for(i in 1:11) {
    for(j in 1:num_splines) {
      chem_df[(1+63*(j-1)):(63*(j)),i] <- as.numeric(chem$fd[[i]]$coefs[j,])
    }
    chem_names <- append(chem_names, paste0("C", i))
  }
  names(chem_df) <- append(chem_names, "category")
  chem_df$category <- ifelse(chem$groupd==1 | chem$groupd==10, 
                             as.character(chem$groupd), "other")
  dev.new()
  if (guide) {
    animate_xy(chem_df[, 1:11], col = rep(chem_df$category, num_splines), 
               fps = 30, tour_path = guided_tour(holes()))
  } else {
    animate_xy(chem_df[, 1:11], col = rep(chem_df$category, num_splines), 
               fps = 30)
  }
}



# Load Data
chem <- merge_chemical_csv(allowed_classes = 1:10, nbasis = 20)



# Visualize the 5th Coefficient
animate_single_coeff(5, F)



# Visualize the 5th and 15th Coefficient
animate_multiple_coeff(c(5, 15), F)



# Visualize all coefficients
animate_all_coeff(20, F)
