#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")
set.seed(999)



#################### IMPORT PACKAGES ###########################
# library("dplyr")
library("tidyverse")
# library("readr")
library("fda")
library("fda.usc")



#################### GLOBAL VARIABLES ##########################
DETECT_OUTLIERS <- FALSE
OUTLIER_TRIM <- 0.1



#################### DEFINE FUNCTIONS ##########################



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
plot_matrix <- data.matrix(plot_df)
# plot_fdata <- fdata.cen(plot_matrix)$Xcen
plot_fdata <- fdata(plot_matrix)
pdf("class_1.pdf")
plot.fdata(
     plot_fdata, 
     type = "l", 
     col = "red",
     main = "Class 1",
     xlab = "Time",
     ylab = "Value"
)

# Plot class 2
drops <- c("X1")
plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
# plot_matrix <- t(data.matrix(plot_df))
plot_matrix <- data.matrix(plot_df)
# plot_fdata <- fdata.cen(plot_matrix)$Xcen
plot_fdata <- fdata(plot_matrix)
pdf("class_2.pdf")
plot.fdata(
     plot_fdata, 
     type = "l", 
     col = "black",
     main = "Class 2",
     xlab = "Time",
     ylab = "Value"
)

# Plot both classes
drops <- c("X1")
plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
plot_matrix <- data.matrix(plot_df)
plot_fdata <- fdata(plot_matrix)
pdf("all_classes.pdf")
plot.fdata(
     plot_fdata, 
     type = "l", 
     col = "red",
     main = "All classes",
     xlab = "Time",
     ylab = "Value"
)
drops <- c("X1")
plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
# plot_matrix <- t(data.matrix(plot_df))
plot_matrix <- data.matrix(plot_df)
plot_fdata <- fdata(plot_matrix)
lines(plot_fdata, col="black")

if (DETECT_OUTLIERS) {

    # Search for outliers in class 1
    drops <- c("X1")
    ecg_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    ecg_matrix <- data.matrix(ecg_df)
    ecg_fdata <- fdata(ecg_matrix)
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim=OUTLIER_TRIM,)
    num_of_outliers <- 0
    for (ele in ecg_outliers$outliers) {
        num_of_outliers <- num_of_outliers + 1
    }
    cat("Number of outliers in Class 1:", num_of_outliers)
    pdf("class_1_outliers.pdf")
    plot.fdata(
               ecg_fdata,
               type = "l",
               col = "red",
               main = "Outliers in Class 1",
               xlab = "Time",
               ylab = "Value",
    )
    for (otlr in ecg_outliers$outliers) {
        otlr <- strtoi(otlr)
        lines(ecg_fdata$data[otlr, ], col="blue")
    }

    # Search for outliers in class 2
    drops <- c("X1")
    ecg_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    ecg_matrix <- data.matrix(ecg_df)
    ecg_fdata <- fdata(ecg_matrix)
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim=OUTLIER_TRIM)
    num_of_outliers <- 0
    for (ele in ecg_outliers$outliers) {
        num_of_outliers <- num_of_outliers + 1
    }
    cat("\nNumber of outliers in Class 2:", num_of_outliers)
    pdf("class_2_outliers.pdf")
    plot.fdata(
               ecg_fdata,
               type = "l",
               col = "black",
               main = "Outliers in Class 2",
               xlab = "Time",
               ylab = "Value"
    )
    for (otlr in ecg_outliers$outliers) {
        otlr <- strtoi(otlr)
        lines(ecg_fdata$data[otlr, ], col="blue")
    }
}
