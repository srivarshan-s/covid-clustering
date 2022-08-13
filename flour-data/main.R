#################### ENIRONMENT SETUP ##########################
setwd("~/Documents/flour-data/")
source("tfunHDDC.R")
source("cfunHDDC.R")



#################### IMPORT PACKAGES ###########################

library("tidyverse")
library("fda")
library("fda.usc")
library("funHDDC")
library("readxl")
library("gtools")
library("RColorBrewer") # color pallette library



#################### FUNCTIONS #################################

# Function that sets seed to reduce randomness
set_seed <- function() {
  set.seed(14071)
}


# Function that loads the functional data and return a list with the functional
# data and the labels
# Parameters:
#       type      => pass the type of basis, "fourier" or "bspline"
#       nbasis    => pass the nbasis or nsplines for the functional data
#       interval  => pass a list where ll is the lower limit, ul is the upper
#                    limit and len is the length of the data
functional_data <- function(type, nbasis, interval) {
  data <- read_excel("data/farines_tri_codage_article_leved_132.xls")
  labels <- data$Qualite
  drops <- c("NomFar", "Qualite")
  data <- data[, !(names(data) %in% drops)]
  fdata <- c()
  if (type == "fourier") {
    set_seed()
    basis <- fda::create.fourier.basis(
      rangeval = c(interval$ll, interval$ul),
      nbasis = nbasis
    )
  } else if (type == "bspline") {
    set_seed()
    basis <- create.bspline.basis(
      rangeval = c(interval$ll, interval$ul),
      nbasis = nbasis
    )
  } else {
    stop("Invalid basis!")
  }
  fdata <- smooth.basis(
    argvals = seq(interval$ll, interval$ul, length.out = interval$len),
    y = t(data.matrix(data)),
    fdParobj = basis
  )$fd
  return(list(fdata=fdata, labels = labels, data = data))
}


# Function that detects the outliers present in the data and returns a list
# with the functional data, labels and outliers
# Parameters:
#       data      => pass the data object returned by the functional_data()
#                    function
#       trim      => pass the alpha of trimming
detect_outliers <- function(data, trim) {
  outliers <- rep(0, ncol(data$fdata$coefs))
  outlier_labels <- c()
  for (ele in unique(data$labels)) {
    data_subset <- data$fdata[data$labels == ele]
    set_seed()
    data_outliers <- outliers.depth.trim(data_subset, trim = trim)
    cat("Class", ele, "outliers:\n")
    print(data_outliers$outliers)
    outlier_labels <- append(outlier_labels, data_outliers$outliers)
  }
  set_seed()
  data_outliers <- outliers.depth.trim(data$fdata, trim = trim)
  cat("Entire data outliers:\n")
  print(data_outliers$outliers)
  new_labels <- data$labels
  for (otlr in outlier_labels) {
    otlr <- strsplit(otlr, split = "rep")[[1]][2]
    otlr <- strtoi(otlr)
    outliers[otlr] <- 1
  }
  data <- list(fdata=data$fdata, labels=data$labels, outliers=outliers, 
               data = data$data)
  return(data)
}


# Function that returns the CCR when given the actual labels and the predicted
# labels
# Parameters:
#       labels    => array containing the actual labels
#       pred      => array containing the predicted labels
find_ccr <- function(labels, pred) {
  cf_matrix <- table(labels, pred)
  ccr <- sum(diag(cf_matrix))
  ccr <- ccr / sum(cf_matrix)
  return(ccr)
}


# Function that finds the optimal numbering for the cluster when given actual
# and predicted labels as parameters
# Parameters:
#       labels    => array containing the actual labels
#       pred      => array containing the predicted labels
find_optimal_labels <- function(labels, pred) {
  max_ccr <- 0
  new_vals <- unique(labels)
  num_class <- length(unique(new_vals))
  max_vals <- new_vals
  for (i in 1:factorial(num_class)) {
    vals <- permutations(n = num_class, r = num_class, v = new_vals)[i, ]
    pred_temp <- vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    if (ccr > max_ccr) {
      max_vals <- vals
      max_ccr <- ccr
    }
  }
  return(max_vals[as.factor(pred)])
}


# Function that prints the number of misclassified labels and the number of
# outliers in the data that are misclassified
# Parameters:
#       result    => list return by the clustering algorithms
#       data      => data object (list) returned by the detect_outliers() 
#                    function
find_misclassified_labels <- function(result, data) {
  misclassified_labels <- c()
  for (idx in 1:length(data$labels)) {
    if (data$labels[idx] != result$class[idx]) {
      misclassified_labels <- append(misclassified_labels, idx)
    }
  }
  cat("Number of misclassified labels:", length(misclassified_labels), 
      "/",  length(data$labels), "\n")
  misclassified_outliers <- c()
  for (idx in misclassified_labels) {
    if (data$outliers[idx] == 1) {
      misclassified_outliers <- append(misclassified_outliers, idx)
    }
  }
  cat("Number of misclassified outliers:", length(misclassified_outliers),
      "/",  sum(data$outliers), "\n")
}


# Function that prints the eta values returned by the cfunHDDC algorithm
# (this does not work for other algorithms)
# Parameters:
#       result    => list return by the cfunHDDC algorithm
find_eta_values <- function(result) {
  eta_vals <- c()
  for (eta in result$etax) {
    eta_vals <- append(eta_vals, eta)
  }
  cat("Eta values:", eta_vals, "\n")
}


# Function that prints the outliers detected by the cfunHDDC algorithm
# Parameters:
#       labels    => the actual labels
#       outliers  => the outlier parameter from the list returned by the
#                    cfunHDDC algorithm
cfunhddc_outliers <- function(labels, outliers) {
  class_outliers <- rep(0, length(unique(labels)))
  for (idx in 1:length(labels)) {
    if (outliers[idx] == 0) {
      class_outliers[labels[idx]] <- class_outliers[labels[idx]] + 1
    }
  }
  for (idx in 1:length(class_outliers)) {
    cat(
      "Number of outliers in Class", idx, ":", class_outliers[idx],
      "/", length(labels[labels == idx]), "\n"
    )
    
  }
}

# Function that takes in the data object returned by functional_data() or
# detect_outliers() functions and performs clustering based on a few parameters
# that are also passed to the function (using the data object returned by the
# detect_outliers() functions enables the function to additionally print the
# misclassified labels and outliers) and returns the result object
# Parameters:
#       data      => data object (list) returned by the functional_data() or 
#                    detect_outliers() function
#       alg       => clustering algorithm to use, funHDDC, cfunHDDC or tfunHDDC
#       K         => number of clusters
#       init      => array containing the nationalizations, kmeans or random
#       nb.rep    => number of repetitions
#       itermax   => maximum number of iterations
#       threshold => threshold for the clustering algorithm
#       models    => array conatining the models to run
#       dfstart   => parameter for tfunHDDC
#       dfupdate  => parameter for tfunHDDC
#       dconstr   => parameter for tfunHDDC
#       alphamin  => parameter for cfunHDDC
cluster_data <- function(data, alg = "funHDDC", K = 2, init = c("kmeans"), 
                         nb.rep = 1, itermax = 200, threshold = c(0.1),
                         models = c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk",
                                    "AkBQkDk", "ABkQkDk", "ABQkDk"),
                         dfstart = 50, dfupdate = c("approx"), 
                         dconstr = c("no"), alphamin = c(0.8), rand_seed = T) {
  result <- NULL
  if (alg == "funHDDC") {
    print("Running funHDDC")
    for (.init in init) {
      for (.threshold in threshold) {
        cat("\n\n\n")
        if (rand_seed) {set_seed()}
        result <- funHDDC(
          data$fdata,
          K = K,
          init = .init,
          threshold = .threshold,
          model = models,
          itermax = itermax,
          nb.rep = nb.rep
        )
        if (K == length(unique(data$labels))) {
          result$class <- find_optimal_labels(data$labels, result$class)
          cf_matrix <- table(data$labels, result$class)
          ccr <- find_ccr(data$labels, result$class)
          print(cf_matrix)
          cat(
            "threshold", .threshold,
            "init:", .init, "ccr:", ccr, "\n"
          )
          if (!is.null(data$outliers)) {
            find_misclassified_labels(result, data)
          }
        } else {
          print("Number of clusters is not the same of number of classes!")
        }
      }
    }
  } else if (alg == "tfunHDDC") {
    print("Running tfunHDDC")
    for (.init in init) {
      for (.threshold in threshold) {
        for (.dfupdate in dfupdate) {
          for (.dconstr in dconstr) {
            cat("\n\n\n")
            if (rand_seed) {set_seed()}
            result <- tfunHDDC(
              data$fdata,
              K = K,
              init = .init,
              threshold = .threshold,
              dfstart = dfstart,
              dfupdate = .dfupdate,
              dconstr = .dconstr,
              model = models,
              itermax = itermax,
              nb.rep = nb.rep
            )
            if (K == length(unique(data$labels))) {
              result$class <- find_optimal_labels(data$labels, result$class)
              cf_matrix <- table(data$labels, result$class)
              ccr <- find_ccr(data$labels, result$class)
              print(cf_matrix)
              cat(
                "threshold", .threshold,
                "init:", .init,
                "dfupdate:", .dfupdate,
                "dconstr:", .dconstr,
                "ccr:", ccr, "\n"
              )
              if (!is.null(data$outliers)) {
                find_misclassified_labels(result, data)
              }
            } else {
              print("Number of clusters is not the same of number of classes!")
            }
          }
        }
      }
    }
  } else if (alg == "cfunHDDC") {
    print("Running cfunHDDC")
    for (.init in init) {
      for (.threshold in threshold) {
        for (.alphamin in alphamin) {
          cat("\n\n\n")
          if (rand_seed) {set_seed()}
          result <- cfunHDDC(
            data$fdata,
            K = K,
            init = .init,
            threshold = .threshold,
            alphamin = .alphamin,
            model = models,
            itermax = itermax,
            nb.rep = nb.rep
          )
          if (K == length(unique(data$labels))) {
            result$class <- find_optimal_labels(data$labels, result$class)
            cf_matrix <- table(data$labels, result$class)
            ccr <- find_ccr(data$labels, result$class)
            print(cf_matrix)
            cat(
              "threshold", .threshold,
              "init:", .init,
              "alphamin:", .alphamin,
              "ccr:", ccr, "\n"
            )
            if (!is.null(data$outliers)) {
              find_misclassified_labels(result, data)
            }
            find_eta_values(result)
            cfunhddc_outliers(data$labels, result$outlier)
          } else {
            print("Number of clusters is not the same of number of classes!")
          }
        }
      }
    }
  } else {
    stop("Invalid algorithm!")
  }
  return(result)
}


CL <- brewer.pal(n = 8, name = "Set1") # color to use for graphs


# Function that plots the original data
# Parameters:
#       data      => object returned by functional_data() or detect_outliers()
plot_data <- function(data) {
  par(mfrow = c(2,2))
  matplot(t(data$data), col = CL[data$labels], main = "a",
          type = "l", lty = 1, xlab = "Time")
  data1 <- data$data[which(data$labels == 1), ]
  matplot(t(as.matrix(data1)), col = CL[1], main = "b",
       type = "l", lty = 1, xlab = "Time")
  data2 <- data$data[which(data$labels == 2), ]
  matplot(t(as.matrix(data2)), col = CL[2], main = "c",
       type = "l", lty = 1, xlab = "Time")
  data3 <- data$data[which(data$labels == 3), ]
  matplot(t(as.matrix(data3)), col = CL[3], main = "d",
       type = "l", lty = 1, xlab = "Time")
}


# Function that plots the functional data
# Parameters:
#       data      => object returned by functional_data() or detect_outliers()
plot_func_data <- function(data) {
  par(mfrow = c(2,2))
  plot(data$fdata, col = CL[data$labels], main = "a", xlab = "Time")
  data1 <- data$fdata[which(data$labels == 1), ]
  plot(data1, col = CL[1], main = "b", xlab = "Time")
  data2 <- data$fdata[which(data$labels == 2), ]
  plot(data2, col = CL[2], main = "c", xlab = "Time")
  data3 <- data$fdata[which(data$labels == 3), ]
  plot(data3, col = CL[3], main = "d", xlab = "Time")
}


# Function that plots the original data with the outliers
# Parameters:
#       data      => object returned by detect_outliers()
plot_contm_data <- function(data) {
  if (is.null(data$outliers)) {stop("Run outlier detection first!")}
  par(mfrow = c(2,2))
  new_labels <- data$labels
  for (idx in 1:length(new_labels)) {
    if (data$outliers[idx] == 1) {
      new_labels[idx] <- 4
    }
  }
  matplot(t(data$data), col = CL[new_labels], main = "a",
          type = "l", lty = 1, xlab = "Time")
  data1 <- data$data[which(data$labels == 1), ]
  new_labels1 <- new_labels[which(data$labels==1)]
  matplot(t(as.matrix(data1)), col = CL[new_labels1], main = "b",
          type = "l", lty = 1, xlab = "Time")
  data2 <- data$data[which(data$labels == 2), ]
  new_labels2 <- new_labels[which(data$labels==2)]
  matplot(t(as.matrix(data2)), col = CL[new_labels2], main = "c",
          type = "l", lty = 1, xlab = "Time")
  data3 <- data$data[which(data$labels == 3), ]
  new_labels3 <- new_labels[which(data$labels==3)]
  matplot(t(as.matrix(data3)), col = CL[new_labels3], main = "d",
          type = "l", lty = 1, xlab = "Time")
}


# Function that plots the functional data with the outliers
# Parameters:
#       data      => object returned by detect_outliers()
plot_contm_func_data <- function(data) {
  if (is.null(data$outliers)) {stop("Run outlier detection first!")}
  par(mfrow = c(2,2))
  new_labels <- data$labels
  for (idx in 1:length(new_labels)) {
    if (data$outliers[idx] == 1) {
      new_labels[idx] <- 4
    }
  }
  plot(data$fdata, col = CL[new_labels], main = "a", xlab = "Time")
  data1 <- data$fdata[which(data$labels == 1), ]
  new_labels1 <- new_labels[which(data$labels==1)]
  plot(data1, col = CL[new_labels1], main = "b", xlab = "Time")
  data2 <- data$fdata[which(data$labels == 2), ]
  new_labels2 <- new_labels[which(data$labels==2)]
  plot(data2, col = CL[new_labels2], main = "c", xlab = "Time")
  data3 <- data$fdata[which(data$labels == 3), ]
  new_labels3 <- new_labels[which(data$labels==3)]
  plot(data3, col = CL[new_labels3], main = "d", xlab = "Time")
}


# Function that plots the clusters with the outliers
# Parameters:
#       result    => object returned by cluster_data()
#       data      => object returned by detect_outliers()
plot_cluster <- function(result, data) {
  if (is.null(data$outliers)) {stop("Run outlier detection first!")}
  par(mfrow = c(2,2))
  new_labels <- result$class
  for (idx in 1:length(new_labels)) {
    if (data$outliers[idx] == 1) {
      new_labels[idx] <- 4
    }
  }
  matplot(t(data$data), col = CL[result$class], main = "a",
          type = "l", lty = 1, xlab = "Time")
  data1 <- data$data[which(result$class == 1), ]
  new_labels1 <- new_labels[which(result$class==1)]
  matplot(t(as.matrix(data1)), col = CL[new_labels1], main = "b",
          type = "l", lty = 1, xlab = "Time")
  data2 <- data$data[which(result$class == 2), ]
  new_labels2 <- new_labels[which(result$class==2)]
  matplot(t(as.matrix(data2)), col = CL[new_labels2], main = "c",
          type = "l", lty = 1, xlab = "Time")
  data3 <- data$data[which(result$class == 3), ]
  new_labels3 <- new_labels[which(result$class==3)]
  matplot(t(as.matrix(data3)), col = CL[new_labels3], main = "d",
          type = "l", lty = 1, xlab = "Time")
}


plot_func_cluster <- function(result, data) {
  if (is.null(data$outliers)) {stop("Run outlier detection first!")}
  par(mfrow = c(2,2))
  new_labels <- result$class
  for (idx in 1:length(new_labels)) {
    if (data$outliers[idx] == 1) {
      new_labels[idx] <- 4
    }
  }
  plot(data$fdata, col = CL[result$class], main = "a", xlab = "Time")
  data1 <- data$fdata[which(result$class == 1), ]
  new_labels1 <- new_labels[which(result$class==1)]
  plot(data1, col = CL[new_labels1], main = "b", xlab = "Time")
  data2 <- data$fdata[which(result$class == 2), ]
  new_labels2 <- new_labels[which(result$class==2)]
  plot(data2, col = CL[new_labels2], main = "c", xlab = "Time")
  data3 <- data$fdata[which(result$class == 3), ]
  new_labels3 <- new_labels[which(result$class==3)]
  plot(data3, col = CL[new_labels3], main = "d", xlab = "Time")
}



#################### MAIN CODE #################################

data <- functional_data(type = "fourier", nbasis = 11, interval = list(
  ll = 0, ul = 1, len = 241))


data <- detect_outliers(data = data, trim = 0.1)


# result <- cluster_data(data = data, alg = "funHDDC", init = "kmeans", K = 3,
#              threshold = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4))


# result <- cluster_data(data = data, alg = "cfunHDDC", init = "kmeans", K = 3,
#              threshold = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4),
#              alphamin = c(0.7, 0.8, 0.85, 0.9, 0.95))


# result <- cluster_data(data = data, alg = "tfunHDDC", init = "kmeans", K = 3,
#              threshold = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4),
#              dfupdate = c("numeric", "approx"),
#              dconstr = c("yes", "no"))

# result <- cluster_data(data = data, alg = "funHDDC", init = "kmeans", K = 3, nb.rep = 50,
#              threshold = c(0.3, 0.4), rand_seed = F)

# result <- cluster_data(data = data, alg = "tfunHDDC", init = "kmeans", K = 3, nb.rep = 10,
#              threshold = c(0.3, 0.4),
#              dfupdate = c("numeric", "approx"),
#              dconstr = c("yes", "no"))

# result <- cluster_data(data = data, alg = "cfunHDDC", init = "kmeans", K = 3, nb.rep = 10,
#              threshold = c(0.1, 0.2, 0.3, 0.4),
#              alphamin = c(0.7, 0.8, 0.85, 0.9, 0.95))
