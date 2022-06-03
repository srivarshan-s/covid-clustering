#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")
source("~/Documents/functional-data-clustering/src/tfunHDDC.R")



#################### IMPORT PACKAGES ###########################
library("tidyverse")
library("fda")
library("fda.usc")
library("funHDDC")



#################### GLOBAL VARIABLES ##########################
DETECT_OUTLIERS <- FALSE
OUTLIER_TRIM <- 0.1
FOURIER_BASIS <- FALSE # TRUE -> fourier FALSE -> bspline 
NBASIS_FOURIER <- 11
NSPLINE_BSPLINE <- 20
ITER_MAX <- 200
MODELS <- c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "AkBQkDk",
            "ABkQkDk", "ABQkDk")



#################### DEFINE FUNCTIONS ##########################

set_seed <- function() {
    set.seed(14071)
}

functional_data <- function(data) {
    if (FOURIER_BASIS) {
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
    } else {
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
}

extract_labels <- function(df) {
    labels <- data.matrix(df[, c("X1")])
    for (idx in 1:length(labels)) {
        if (labels[idx] == -1) {
            # Mapping  1 -> 1
            #         -1 -> 2
            labels[idx] <- 2
        }
    }
    return(labels)
}

change_labels <- function(labels) {
    for (idx in 1:length(labels)) {
        if (labels[idx] == 1) {
            labels[idx] <- 2
        } else {
            labels[idx] <- 1
        }
    }
    return(labels)
}



#################### MAIN CODE #################################

# Read train data
df_train <- read_tsv("data/ECG200_TRAIN.tsv", col_names = FALSE)
print("Training data")
print(df_train)

# Read test data
df_test <- read_tsv("data/ECG200_TEST.tsv", col_names = FALSE)
print("Testing data")
print(df_test)

# Merging the two dataframes
df <- rbind(df_train, df_test)
print("Merged data")
print(df)

# Extract the labels
labels <- extract_labels(df)

# Extracting the two classes
df_class_1 <- filter(df, X1 == "-1")
df_class_2 <- filter(df, X1 == "1")
print("Class 1")
print(df_class_1)
print("Class 2")
print(df_class_2)

# Plot class 1
drops <- c("X1")
plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
plot_fdata <- functional_data(plot_df)
pdf("class_1.pdf")
plot.fd(
        plot_fdata, 
        col = "red",
        main = "Class 1",
        xlab = "Time",
        ylab = "Value"
)

# Plot class 2
drops <- c("X1")
plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
plot_fdata <- functional_data(plot_df)
pdf("class_2.pdf")
plot.fd(
        plot_fdata,
        col = "black",
        main = "Class 2",
        xlab = "Time",
        ylab = "Value"
)

# Plot both classes
drops <- c("X1")
plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
plot_fdata <- functional_data(plot_df)
pdf("all_classes.pdf")
plot.fd(
        plot_fdata, 
        col = "red",
        main = "All classes",
        xlab = "Time",
        ylab = "Value"
)
drops <- c("X1")
plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
plot_fdata <- functional_data(plot_df)
lines(plot_fdata, col = "black")

# Outlier detection
if (DETECT_OUTLIERS) {

    # Search for outliers in class 1
    print("Detecting outliers in Class 1.....")
    drops <- c("X1")
    ecg_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim = OUTLIER_TRIM)
    num_of_outliers <- 0
    cat("Number of outliers in Class 1:",
        length(ecg_outliers$outliers), "\n")
    pdf("class_1_outliers.pdf")
    plot.fd(
            ecg_fdata,
            col = "red",
            main = "Outliers in Class 1",
            xlab = "Time",
            ylab = "Value"
    )
    for (otlr in ecg_outliers$outliers) {
        lines(
              ecg_fdata[ecg_fdata$fdnames$reps == otlr],
              col = "blue"
        )
    }

    # Search for outliers in class 2
    print("Detecting outliers in Class 2.....")
    drops <- c("X1")
    ecg_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    set_seed()
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim = OUTLIER_TRIM)
    num_of_outliers <- 0
    cat("Number of outliers in Class 2:",
        length(ecg_outliers$outliers), "\n")
    pdf("class_2_outliers.pdf")
    plot.fd(
            ecg_fdata,
            col = "black",
            main = "Outliers in Class 2",
            xlab = "Time",
            ylab = "Value"
    )
    for (otlr in ecg_outliers$outliers) {
        lines(
              ecg_fdata[ecg_fdata$fdnames$reps == otlr],
              col = "blue"
        )
    }
}

# # funHDDC algorithm
# print("Running funHDDC algorithm.....")
# drops <- c("X1")
# ecg_df <- df[, !(names(df) %in% drops)]
# ecg_fdata <- functional_data(ecg_df)
# set_seed()
# result <- funHDDC(
#                   ecg_fdata,
#                   K = 2,
#                   init = "random", # 'random', 'kmeans'
#                   threshold = 0.001,
#                   model = MODELS,
#                   itermax = ITER_MAX,
#                   nb.rep = 50
# )
# if (FOURIER_BASIS) {
#     pdf("funHDDC_fourier_clusters.pdf")
# } else {
#     pdf("funHDDC_bspline_clusters.pdf")
# }
# plot.fd(ecg_fdata, col = result$class, lwd = 2, lty = 1)
# cf_matrix <- table(labels, result$class)
# if (ccr < 1 - ccr) {
#     labels <- change_labels(labels)
#     cf_matrix <- table(labels, result$class)
#     ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# }
# print(cf_matrix)
# ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# cat("The correct classification rate:", ccr * 100, "%\n")

# funHDDC gridsearch
print("Running funHDDC gridsearch.....")
drops <- c("X1")
ecg_df <- df[, !(names(df) %in% drops)]
ecg_fdata <- functional_data(ecg_df)
GRIDSEARCH_INITS <- c("kmeans", "random")
GRIDSEARCH_THRESHOLDS <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
BEST_CCR <- 0
BEST_INIT <- ""
BEST_THRESHOLD <- 0
for (init in GRIDSEARCH_INITS) {
    for (threshold in GRIDSEARCH_THRESHOLDS) {
        cat("\n\n\n")
        set_seed()
        result <- funHDDC(
                          ecg_fdata,
                          K = 2,
                          init = init,
                          threshold = threshold,
                          model = MODELS,
                          itermax = ITER_MAX,
        )
        cf_matrix <- table(labels, result$class)
        ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
        if (ccr < 1 - ccr) {
            labels <- change_labels(labels)
            cf_matrix <- table(labels, result$class)
            ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
        }
        print(cf_matrix)
        cat("threshold", threshold, 
            "init:", init, "ccr:", ccr, "\n")
        if (ccr >= BEST_CCR) {
            BEST_CCR <- ccr
            BEST_INIT <- init
            BEST_THRESHOLD <- threshold
        }
    }
}
cat("Best Init:", BEST_INIT, "\n")
cat("Best Threshold:", BEST_THRESHOLD, "\n")
cat("Highest CCR:", BEST_CCR, "\n")

# Find the number of misclassified data
misclassified_labels <- c()
for (idx in 1:length(labels)) {
    if (labels[idx] != result$class[idx]) {
        misclassified_labels <- append(
                                       misclassified_labels, 
                                       paste("rep", as.character(idx), sep = "")
        )
    }
}
cat("Number of misclassified labels:", length(misclassified_labels), "\n")

# Find the number of outliers that are misclassified

# # tfunHDDC algorithm
# print("Running tfunHDDC algorithm.....")
# drops <- c("X1")
# ecg_df <- df[, !(names(df) %in% drops)]
# ecg_fdata <- functional_data(ecg_df)
# set_seed()
# result <- tfunHDDC(
#   ecg_fdata,
#   K = 2,
#   init = "kmeans", # 'random', 'kmeans'
#   threshold = 0.01,
#   model = GRIDSEARCH_MODELS,
#   itermax = FUNHDDC_ITER_MAX,
#   nb.rep = 1,
#   dfstart=50, # Don't go more than 100; similar to funHDDC; keep it at 50
#   dfupdate = "numeric", # "approx", "numeric"
#   dconstr = "no" # "yes", "no"
# )
# cf_matrix <- table(labels, result$class)
# ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# if (ccr < 1 - ccr) {
#   labels <- change_labels(labels)
#   cf_matrix <- table(labels, result$class)
#   ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# }
# print(cf_matrix)
# cat("The correct classification rate:", ccr * 100, "%\n")
