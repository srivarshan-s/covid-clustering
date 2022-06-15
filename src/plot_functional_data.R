#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")



#################### IMPORT PACKAGES ###########################
library("tidyverse")
library("fda")
library("fda.usc")
library("funHDDC")



#################### GLOBAL VARIABLES ##########################
NBASIS_FOURIER <- 41
NSPLINE_BSPLINE <- 38



#################### DEFINE FUNCTIONS ##########################

set_seed <- function() {
    set.seed(14071)
}

fourier_basis <- function(data) {
    set_seed()
    basis <- fda::create.fourier.basis(
                                       rangeval = c(0, 95),
                                       nbasis = NBASIS_FOURIER
    )
    ecg_fdata <- smooth.basis(
                              argvals = seq(0, 95, length.out = 96),
                              y = t(data.matrix(data)),
                              fdParobj = basis
                              )$fd
    return(ecg_fdata)
}

bspline_basis <- function(data) {
    set_seed()
    basis <- create.bspline.basis(
                                  rangeval = c(0, 95), 
                                  nbasis = NSPLINE_BSPLINE
    )
    ecg_fdata <- smooth.basis(
                              argvals = seq(0, 95, length.out = 96),
                              y = t(data.matrix(data)),
                              fdParobj = basis
                              )$fd
    return(ecg_fdata)
}



#################### MAIN CODE #################################

# Read train data
df_train <- read_tsv("data/ECG200_TRAIN.tsv", col_names = FALSE)

# Read test data
df_test <- read_tsv("data/ECG200_TEST.tsv", col_names = FALSE)

# Merging the two dataframes
df <- rbind(df_train, df_test)

# Extract normal data
drops <- c("X1")
plot_df <- df[, !(names(df) %in% drops)]
normal_data <- data.matrix(plot_df)

# Extract fourier data
fourier_data <- fourier_basis(plot_df)

# Extract bspline data
bspline_data <- bspline_basis(plot_df)

# Plot data
for (idx in 151:200) {
    pdf_name <- paste("norm_vs_func_", as.character(idx), ".pdf", sep = "")
    pdf(pdf_name)
    plot.fd(
            fourier_data[idx],
            col = 2,
            xlab = "Time",
            ylab = "Value",
    )
    lines(normal_data[idx, ], col = 1)
    lines(bspline_data[idx], col = 3)
}
