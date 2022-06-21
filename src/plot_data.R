#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")
source("./src/helper_functions.R")
source("./src/tfunHDDC.R")
source("./src/cfunHDDC.R")



#################### IMPORT PACKAGES ###########################
library("tidyverse")
library("fda")
library("fda.usc")
library("funHDDC")



#################### GLOBAL VARIABLES ##########################
NBASIS_FOURIER <- 11
NSPLINE_BSPLINE <- 20
FOURIER_BASIS <- FALSE
PLOT_ORIGINAL_DATA <- FALSE
PLOT_TFUNHDDC_CLUSTERS <- FALSE
PLOT_CFUNHDDC_CLUSTERS <- TRUE
ITER_MAX <- 200
MODELS <- c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "AkBQkDk",
            "ABkQkDk", "ABQkDk")



#################### MAIN CODE #################################

# Read train data
df_train <- read_tsv("data/ECG200_TRAIN.tsv", col_names = FALSE)
# print("Training data")
# print(df_train)

# Read test data
df_test <- read_tsv("data/ECG200_TEST.tsv", col_names = FALSE)
# print("Testing data")
# print(df_test)

# Merging the two dataframes
df <- rbind(df_train, df_test)
# print("Merged data")
# print(df)

# Extract the labels
labels <- extract_labels(df)

# Extracting the two classes
df_class_1 <- filter(df, X1 == "-1")
df_class_2 <- filter(df, X1 == "1")
# print("Class 1")
# print(df_class_1)
# print("Class 2")
# print(df_class_2)

if (PLOT_ORIGINAL_DATA) {

    # Plot class 1 of original data
    drops <- c("X1")
    plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    pdf("original_data_class_1.pdf")
    plot(
         as.numeric(plot_df[1, ]),
         type = "l",
         col = "red",
         main = "Class 1",
         xlab = "Time",
         ylab = "Value"
    )
    for (idx in 2:67) {
        lines(as.numeric(plot_df[idx, ]), col = "red")
    }

    # Plot class 2 of original data
    drops <- c("X1")
    plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    pdf("original_data_class_2.pdf")
    plot(
         as.numeric(plot_df[1, ]),
         type = "l",
         col = "black",
         main = "Class 2",
         xlab = "Time",
         ylab = "Value"
    )
    for (idx in 2:133) {
        lines(as.numeric(plot_df[idx, ]), col = "black")
    }

    # Plot both classes of original data
    drops <- c("X1")
    plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    pdf("original_data_all_classes.pdf")
    plot(
         as.numeric(plot_df[1, ]),
         type = "l",
         col = "red",
         main = "All classes",
         xlab = "Time",
         ylab = "Value"
    )
    for (idx in 2:67) {
        lines(as.numeric(plot_df[idx, ]), col = "red")
    }
    drops <- c("X1")
    plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    for (idx in 1:133) {
        lines(as.numeric(plot_df[idx, ]), col = "black")
    }
}

if (PLOT_TFUNHDDC_CLUSTERS) {

    drops <- c("X1")
    ecg_df <- df[, !(names(df) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    result <- tfunHDDC(
                       ecg_fdata,
                       K = 2,
                       init = "random", # 'random', 'kmeans'
                       threshold = 0.001,
                       model = MODELS,
                       itermax = ITER_MAX,
                       nb.rep = 20,
                       dfstart=50, # Don't go more than 100; similar to funHDDC; keep it at 50
                       dfupdate = "numeric", # "approx", "numeric"
                       dconstr = "yes" # "yes", "no"
    )

    labels <- extract_labels(df)
    cf_matrix <- table(labels, result$class)
    ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
    if (ccr < 1 - ccr) {
        labels <- change_labels(labels)
        cf_matrix <- table(labels, result$class)
        ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
    }
    print(ccr)
    pdf("tfunHDDC_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class)
    pdf("tfunHDDC_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1)
    pdf("tfunHDDC_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2)
}

if (PLOT_CFUNHDDC_CLUSTERS) {

    drops <- c("X1")
    ecg_df <- df[, !(names(df) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    result <- cfunHDDC(
                       ecg_fdata,
                       K = 2,
                       init = "kmeans", # 'random', 'kmeans'
                       threshold = 0.1,
                       model = MODELS,
                       itermax = ITER_MAX,
                       nb.rep = 20,
                       alphamin = 0.85 # Ideally between 0.8 and 0.95
    )
    labels <- extract_labels(df)
    cf_matrix <- table(labels, result$class)
    ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
    if (ccr < 1 - ccr) {
        labels <- change_labels(labels)
        cf_matrix <- table(labels, result$class)
        ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
    }
    print(ccr)
    pdf("cfunHDDC_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class)
    pdf("cfunHDDC_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1)
    pdf("cfunHDDC_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2)
}
