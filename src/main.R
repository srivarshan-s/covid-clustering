#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")
set.seed(999)



#################### IMPORT PACKAGES ###########################
library("tidyverse")
library("fda")
library("fda.usc")



#################### GLOBAL VARIABLES ##########################
DETECT_OUTLIERS <- FALSE
OUTLIER_TRIM <- 0.1
FOURIER_BASIS <- TRUE # Do not change to FALSE, as code breaks
NSPLINES <- 9



#################### DEFINE FUNCTIONS ##########################

functional_data <- function(data, nsplines=NSPLINES) {
    if (FOURIER_BASIS) {
        basis <- create.fourier.basis(rangeval=c(0,95), nbasis=nsplines)
        ecg_fdata <- smooth.basis(
                                  argvals=seq(0,95,length.out=96),
                                  y=t(data.matrix(data)),
                                  fdParobj=basis
        )$fd
        return(ecg_fdata)
    } else {
        return(fdata(data.matrix(data)))
    }
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
lines(plot_fdata, col="black")

# Outlier detection
if (DETECT_OUTLIERS) {

    # Search for outliers in class 1
    print("Detecting outliers in Class 1.....")
    drops <- c("X1")
    ecg_df <- df_class_1[, !(names(df_class_1) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim=OUTLIER_TRIM)
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
              col="blue"
        )
    }

    # Search for outliers in class 2
    print("Detecting outliers in Class 2.....")
    drops <- c("X1")
    ecg_df <- df_class_2[, !(names(df_class_2) %in% drops)]
    ecg_fdata <- functional_data(ecg_df)
    ecg_outliers <- outliers.depth.trim(ecg_fdata, trim=OUTLIER_TRIM)
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
              col="blue"
        )
    }
}
