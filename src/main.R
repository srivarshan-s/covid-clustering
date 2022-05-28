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
n_obs <- 67
time_span <- 100
time <- sort(runif(n_obs, 0, time_span))



#################### DEFINE FUNCTIONS ##########################

# Function to smooth the line using the fda package's smooth.basis
smooth_line <- function(val) {
   # n_obs <- length(val)
   # time_span <- 100
   # time <- sort(runif(n_obs, 0, time_span))
   val_obs <- val + rnorm(n_obs, 0, 0.05)
   times_basis <- seq(0, time_span, 1)
   knots <- c(seq(0, time_span, 5))
   n_knots <- length(knots)
   n_order <- 4
   n_basis <- length(knots) + n_order - 2
   basis <- create.bspline.basis(
      c(min(times_basis), max(times_basis)),
      n_basis,
      n_order,
      knots
   )
   val_obj <- smooth.basis(argvals = time, y = val_obs, fdParobj = basis)
   return(val_obj)
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
plot_matrix <- t(data.matrix(plot_df))
pdf("class_1.pdf", width = 20, height = 8)
plot(
   # time,
   # smooth_line(plot_matrix[1, ]),
   plot_matrix[1, ],
   col = "white",
   type = "l",
   main = "Class 1",
   xlab = "Observations",
   ylab = "Electrical Activity"
)
for (i in 1:96) {
   # lines(plot_matrix[i, ], col="red")
   lines(smooth_line((plot_matrix[i, ])), lwd = 1, col = "blue")
}
