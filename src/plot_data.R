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
FOURIER_BASIS <- TRUE
PLOT_ORIGINAL_DATA <- TRUE
PLOT_TFUNHDDC_FOURIER_CLUSTERS <- FALSE
PLOT_CFUNHDDC_FOURIER_CLUSTERS <- FALSE
PLOT_TFUNHDDC_BSPLINE_CLUSTERS <- FALSE
PLOT_CFUNHDDC_BSPLINE_CLUSTERS <- FALSE
OUTLIER_TRIM <- 0.1
ITER_MAX <- 200
LINE_TYPE <- 1
MODELS <- c(
            "AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "AkBQkDk",
            "ABkQkDk", "ABQkDk"
)



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

    outlier_labels <- c()

    # Search for outliers in class 1
    drops <- c("X1")
    ecg_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim = OUTLIER_TRIM)
    num_of_outliers <- 0
    outlier_labels <- append(outlier_labels, ecg_outliers$outliers)

    # Plot class 1 of original data
    drops <- c("X1")
    plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    pdf("original_data_class_1.pdf")
    A <- c()
    for (ele in ecg_fdata$fdnames$reps) {
        flag <- TRUE
        for (otlr in ecg_outliers$outliers) {
            if (otlr == ele) {
                flag <- FALSE
                A <- append(A, 3)
            }
        }
        if (flag) {
            A <- append(A, 1)
        }
    }
    matplot(
            t(as.matrix(plot_df)),
            type = "l",
            col = A,
            lty = LINE_TYPE,
            main = "Class 1",
            xlab = "Time",
            ylab = "Value",
    )

    # Search for outliers in class 2
    drops <- c("X1")
    ecg_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim = OUTLIER_TRIM)
    num_of_outliers <- 0
    outlier_labels <- append(outlier_labels, ecg_outliers$outliers)

    # Plot class 2 of original data
    drops <- c("X1")
    plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    pdf("original_data_class_2.pdf")
    A <- c()
    for (ele in ecg_fdata$fdnames$reps) {
        flag <- TRUE
        for (otlr in ecg_outliers$outliers) {
            if (otlr == ele) {
                flag <- FALSE
                A <- append(A, 3)
            }
        }
        if (flag) {
            A <- append(A, 2)
        }
    }
    matplot(
            t(as.matrix(plot_df)),
            type = "l",
            col = A,
            lty = LINE_TYPE,
            main = "Class 2",
            xlab = "Time",
            ylab = "Value",
    )

    # Plot both classes of original data
    drops <- c("X1")
    plot_df <- df[, !(names(df) %in% drops)]
    pdf("original_data_all_classes.pdf")
    new_labels <- change_labels(labels)
    for (otlr in outlier_labels) {
        idx <- strtoi(strsplit(otlr, "rep")[[1]][2])
        new_labels[idx] <- 3
    }
    matplot(
            t(as.matrix(plot_df)),
            type = "l",
            lty = LINE_TYPE,
            col = new_labels,
            main = "All Classes",
            xlab = "Time",
            ylab = "Value",
    )
}

if (PLOT_TFUNHDDC_FOURIER_CLUSTERS) {
    FOURIER_BASIS <- TRUE
    drops <- c("X1")
    ecg_df <- df[, !(names(df) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    result <- tfunHDDC(
                       ecg_fdata,
                       K = 2,
                       init = "random", # 'random', 'kmeans'
                       threshold = 0.01,
                       model = MODELS,
                       itermax = ITER_MAX,
                       nb.rep = 20,
                       dfstart = 50,
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
    pdf("tfunHDDC_fourier_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class,
            main = "All Classes",
            xlab = "Time",
            ylab = "Value",)
    pdf("tfunHDDC_fourier_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1,
            main = "Class 1",
            xlab = "Time",
            ylab = "Value",)
    pdf("tfunHDDC_fourier_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2,
            main = "Class 2",
            xlab = "Time",
            ylab = "Value",)
}

if (PLOT_CFUNHDDC_FOURIER_CLUSTERS) {
    FOURIER_BASIS <- TRUE
    drops <- c("X1")
    ecg_df <- df[, !(names(df) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    result <- cfunHDDC(
                       ecg_fdata,
                       K = 2,
                       init = "random", # 'random', 'kmeans'
                       threshold = 0.01,
                       model = MODELS,
                       itermax = ITER_MAX,
                       nb.rep = 20,
                       alphamin = 0.8 # Ideally between 0.8 and 0.95
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
    pdf("cfunHDDC_fourier_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class,
            main = "All Classes",
            xlab = "Time",
            ylab = "Value",)
    pdf("cfunHDDC_fourier_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1,
            main = "Class 1",
            xlab = "Time",
            ylab = "Value",)
    pdf("cfunHDDC_fourier_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2,
            main = "Class 2",
            xlab = "Time",
            ylab = "Value",)
}

if (PLOT_TFUNHDDC_BSPLINE_CLUSTERS) {
    FOURIER_BASIS <- FALSE
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
                       dfstart = 50,
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
    pdf("tfunHDDC_bspline_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class,
            main = "All Classes",
            xlab = "Time",
            ylab = "Value",)
    pdf("tfunHDDC_bspline_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1,
            main = "Class 1",
            xlab = "Time",
            ylab = "Value",)
    pdf("tfunHDDC_bspline_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2,
            main = "Class 2",
            xlab = "Time",
            ylab = "Value",)
}

if (PLOT_CFUNHDDC_BSPLINE_CLUSTERS) {
    FOURIER_BASIS <- FALSE
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
    pdf("cfunHDDC_bspline_all_clusters.pdf")
    plot.fd(ecg_fdata, col = result$class,
            main = "All Classes",
            xlab = "Time",
            ylab = "Value",)
    pdf("cfunHDDC_bspline_cluster_1.pdf")
    plot.fd(ecg_fdata[result$class == 1], col = 1,
            main = "Class 1",
            xlab = "Time",
            ylab = "Value",)
    pdf("cfunHDDC_bspline_cluster_2.pdf")
    plot.fd(ecg_fdata[result$class == 2], col = 2,
            main = "Class 2",
            xlab = "Time",
            ylab = "Value",)
}
