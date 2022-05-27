#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/functional-data-clustering")



#################### IMPORT PACKAGES ###########################
# library("dplyr")
library("tidyverse")
# library("readr")
library("fda.usc")



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
pdf( "class_1.pdf", width = 20, height = 8 )
plot(
     plot_matrix[1, ],
     # lwd=2,
     col = "red",
     type = "l",
     main = "Class 1",
     xlab = "Observations",
     ylab = "Electrical Activity"
)
for (i in 2:96) {
   # lines(plot_matrix[i, ], col="red")
}
