#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")
source("~/Documents/functional-data-clustering/src/tfunHDDC.R")
source("~/Documents/functional-data-clustering/src/cfunHDDC.R")



#################### IMPORT PACKAGES ###########################
library("tidyverse")
library("fda")
library("fda.usc")
library("funHDDC")



#################### GLOBAL VARIABLES ##########################
DETECT_OUTLIERS <- TRUE
OUTLIER_TRIM <- 0.1
FOURIER_BASIS <- TRUE # TRUE -> fourier FALSE -> bspline 
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

find_misclassified_labels <- function(result, labels, outlier_labels) {
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
    cat("Number of outliers that are misclassified:", 
        length(intersect(outlier_labels, misclassified_labels)), "\n")
}

find_eta_values <- function(result) {
    eta_vals <- c(result$etax[1], result$etax[2])
    cat("Eta values:", eta_vals, "\n")
}



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

    outlier_labels <- c()

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
    outlier_labels <- append(outlier_labels, ecg_outliers$outliers)

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
    outlier_labels <- append(outlier_labels, ecg_outliers$outliers)
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
#                   init = "kmeans", # 'random', 'kmeans'
#                   threshold = 0.1,
#                   model = MODELS,
#                   itermax = ITER_MAX,
#                   nb.rep = 1
# )
# if (FOURIER_BASIS) {
#     pdf("funHDDC_fourier_clusters.pdf")
# } else {
#     pdf("funHDDC_bspline_clusters.pdf")
# }
# plot.fd(ecg_fdata, col = result$class, lwd = 2, lty = 1)
# cf_matrix <- table(labels, result$class)
# ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# if (ccr < 1 - ccr) {
#     labels <- change_labels(labels)
#     cf_matrix <- table(labels, result$class)
#     ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
# }
# print(cf_matrix)
# cat("The correct classification rate:", ccr * 100, "%\n")
# find_misclassified_labels(result, labels, outlier_labels)

# # funHDDC gridsearch
# print("Running funHDDC gridsearch.....")
# drops <- c("X1")
# ecg_df <- df[, !(names(df) %in% drops)]
# ecg_fdata <- functional_data(ecg_df)
# GRIDSEARCH_INITS <- c("kmeans", "random")
# GRIDSEARCH_THRESHOLDS <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
# BEST_CCR <- 0
# BEST_INIT <- ""
# BEST_THRESHOLD <- 0
# for (init in GRIDSEARCH_INITS) {
#     for (threshold in GRIDSEARCH_THRESHOLDS) {
#         cat("\n\n\n")
#         set_seed()
#         result <- funHDDC(
#                           ecg_fdata,
#                           K = 2,
#                           init = init,
#                           threshold = threshold,
#                           model = MODELS,
#                           itermax = ITER_MAX,
#                           nb.rep = 20
#         )
#         cf_matrix <- table(labels, result$class)
#         ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
#         if (ccr < 1 - ccr) {
#             labels <- change_labels(labels)
#             cf_matrix <- table(labels, result$class)
#             ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
#         }
#         print(cf_matrix)
#         cat("threshold", threshold, 
#             "init:", init, "ccr:", ccr, "\n")
#         find_misclassified_labels(result, labels, outlier_labels)
#         if (ccr >= BEST_CCR) {
#             BEST_CCR <- ccr
#             BEST_INIT <- init
#             BEST_THRESHOLD <- threshold
#         }
#     }
# }
# cat("\n\n\n")
# cat("Best Init:", BEST_INIT, "\n")
# cat("Best Threshold:", BEST_THRESHOLD, "\n")
# cat("Highest CCR:", BEST_CCR, "\n")

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
#   threshold = 0.001,
#   model = MODELS,
#   itermax = ITER_MAX,
#   nb.rep = 1,
#   dfstart=50, # Don't go more than 100; similar to funHDDC; keep it at 50
#   dfupdate = "numeric", # "approx", "numeric"
#   dconstr = "yes" # "yes", "no"
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
# find_misclassified_labels(result, labels, outlier_labels)

# # tfunHDDC gridsearch
# print("Running tfunHDDC gridsearch.....")
# drops <- c("X1")
# ecg_df <- df[, !(names(df) %in% drops)]
# ecg_fdata <- functional_data(ecg_df)
# GRIDSEARCH_INITS <- c("kmeans", "random")
# GRIDSEARCH_THRESHOLDS <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
# GRIDSEARCH_DFUPDATE <- c("numeric", "approx")
# GRIDSEARCH_DCONSTR <- c("yes", "no")
# BEST_CCR <- 0
# BEST_INIT <- ""
# BEST_THRESHOLD <- 0
# BEST_DFUPDATE <- ""
# BEST_DCONSTR <- ""
# for (dconstr in GRIDSEARCH_DCONSTR) { for (dfupdate in GRIDSEARCH_DFUPDATE) {
#     for (init in GRIDSEARCH_INITS) { for (threshold in GRIDSEARCH_THRESHOLDS) {
#         cat("\n\n\n")
#         set_seed()
#         result <- tfunHDDC(
#                            ecg_fdata,
#                            K = 2,
#                            init = init,
#                            threshold = threshold,
#                            model = MODELS,
#                            itermax = ITER_MAX,
#                            nb.rep = 20,
#                            dfstart = 50,
#                            dfupdate = dfupdate,
#                            dconstr = dconstr
#         )
#         cf_matrix <- table(labels, result$class)
#         ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
#         if (ccr < 1 - ccr) {
#             labels <- change_labels(labels)
#             cf_matrix <- table(labels, result$class)
#             ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2]) / sum(cf_matrix)
#         }
#         print(cf_matrix)
#         cat("threshold", threshold, 
#             "init:", init, 
#             "dfupdate:", dfupdate, 
#             "dconstr:", dconstr,
#             "ccr:", ccr, "\n")
#         find_misclassified_labels(result, labels, outlier_labels)
#         if (ccr >= BEST_CCR) {
#             BEST_CCR <- ccr
#             BEST_INIT <- init
#             BEST_THRESHOLD <- threshold
#             BEST_DFUPDATE <- dfupdate
#             BEST_DCONSTR <- dconstr
#         }
#     } }
# } }
# cat("\n\n\n")
# cat("Best Init:", BEST_INIT, "\n")
# cat("Best Threshold:", BEST_THRESHOLD, "\n")
# cat("Best Update:", BEST_DFUPDATE, "\n")
# cat("Best Constraint:", BEST_DCONSTR, "\n")
# cat("Highest CCR:", BEST_CCR, "\n")

# # cfunHDDC algorithm
# print("Running cfunHDDC algorithm.....")
# drops <- c("X1")
# ecg_df <- df[, !(names(df) %in% drops)]
# ecg_fdata <- functional_data(ecg_df)
# set_seed()
# result <- cfunHDDC(
#   ecg_fdata,
#   K = 2,
#   init = "kmeans", # 'random', 'kmeans'
#   threshold = 0.1,
#   model = MODELS,
#   itermax = ITER_MAX,
#   nb.rep = 1,
#   alphamin = 0.5 # Ideally between 0.8 and 0.95
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
# find_misclassified_labels(result, labels, outlier_labels)
# find_eta_values(result)

# cfunHDDC gridsearch
print("Running cfunHDDC gridsearch.....")
drops <- c("X1")
ecg_df <- df[, !(names(df) %in% drops)]
ecg_fdata <- functional_data(ecg_df)
GRIDSEARCH_INITS <- c("kmeans", "random")
GRIDSEARCH_THRESHOLDS <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
GRIDSEARCH_ALPHAS <- c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
BEST_CCR <- 0
BEST_INIT <- ""
BEST_THRESHOLD <- 0
BEST_ALPHAMIN <- 0
for (init in GRIDSEARCH_INITS) { for (threshold in GRIDSEARCH_THRESHOLDS) {
    for (alphamin in GRIDSEARCH_ALPHAS) {
        cat("\n\n\n")
        set_seed()
        result <- cfunHDDC(
                           ecg_fdata,
                           K = 2,
                           init = init,
                           threshold = threshold,
                           model = MODELS,
                           itermax = ITER_MAX,
                           nb.rep = 20,
                           alphamin = alphamin
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
            "init:", init, 
            "alphamin:", alphamin,
            "ccr:", ccr, "\n")
        find_misclassified_labels(result, labels, outlier_labels)
        find_eta_values()
        if (ccr >= BEST_CCR) {
            BEST_CCR <- ccr
            BEST_INIT <- init
            BEST_THRESHOLD <- threshold
            BEST_ALPHAMIN <- alphamin
        }
    }
} }
cat("\n\n\n")
cat("Best Init:", BEST_INIT, "\n")
cat("Best Threshold:", BEST_THRESHOLD, "\n")
cat("Best Alphamin:", BEST_ALPHAMIN, "\n")
cat("Highest CCR:", BEST_CCR, "\n")
