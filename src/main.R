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

# # Plot class 1
# drops <- c("X1")
# plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
# plot_matrix <- t(data.matrix(plot_df))
# pdf("class_1.pdf", width = 20, height = 8)
# plot(
#    # time,
#    # smooth_line(plot_matrix[1, ]),
#    plot_matrix[1, ],
#    col = "white",
#    type = "l",
#    main = "Class 1",
#    xlab = "Observations",
#    ylab = "Electrical Activity"
# )
# for (i in 1:96) {
#    # lines(plot_matrix[i, ], col="red")
#    lines(smooth_line((plot_matrix[i, ])), lwd = 1, col = "blue")
# }

# # Plot class 2
# drops <- c("X1")
# plot_df <- df_class_2[, !(names(df_class_2) %in% drops)]
# plot_matrix <- t(data.matrix(plot_df))
# print(length(plot_matrix[1, ]))
# pdf("class_2.pdf", width = 20, height = 8)
# plot(
#    # time,
#    # smooth_line(plot_matrix[1, ]),
#    plot_matrix[1, ],
#    col = "white",
#    type = "l",
#    main = "Class 2",
#    xlab = "Observations",
#    ylab = "Electrical Activity"
# )
# for (i in 1:96) {
#    # lines(plot_matrix[i, ], col="red")
#    lines(smooth_line((plot_matrix[i, ])), lwd = 1, col = "blue")
# }

# Plot class 1
drops <- c("X1")
plot_df <- df_class_1[, !(names(df_class_1) %in% drops)]
plot_matrix <- data.matrix(plot_df)
plot_fdata.c <- fdata.cen(plot_matrix)$Xcen
pdf("class_1.pdf")
plot.fdata(
     plot_fdata.c, 
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
plot_fdata.c <- fdata.cen(plot_matrix)$Xcen
pdf("class_2.pdf")
plot.fdata(
     plot_fdata.c, 
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
plot_fdata.c <- fdata.cen(plot_matrix)$Xcen
pdf("all_classes.pdf")
plot.fdata(
     plot_fdata.c, 
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
plot_fdata.c <- fdata.cen(plot_matrix)$Xcen
lines(
     plot_fdata.c, 
     col = "black",
)

