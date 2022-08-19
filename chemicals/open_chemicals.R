library(funHDDC)
library(fda)
library(tidyverse)

merge_chemical_csv <- function(csv_range = 1:65, dim_use = 1:11, max_time = NULL,
                               allowed_classes = c(1, 10), nbasis = 10) {
  setwd("./Raw_CSV")

  chemical_labels <- read_csv("chemical_index.csv")
  # set up initial empty dimensional data
  dim_data <- list()
  keep_ind <- which(chemical_labels$Class %in% allowed_classes)
  csv_range <- csv_range[keep_ind]

  # get the max dimensions
  for (load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name)

    if (is.null(max_time) || max(load_df[, "Time"]) < max_time) {
      max_time <- as.numeric(max(load_df$Time))
    }
  }
  # create a basis function based on the highest time can use
  basis <- create.bspline.basis(rangeval = c(0, max_time), nbasis = nbasis)

  fd_chem <- list()
  for (load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name)
    # get only allowed rows
    load_df <- load_df[load_df$Time <= max_time, ]

    for (dim_c in dim_use) {
      dim_name <- paste0("C", dim_c)
      fd_temp <- smooth.basis(load_df$Time,
        y = load_df[[dim_name]],
        fdParobj = basis
      )$fd
      if (length(fd_chem) < dim_c) {
        fd_chem[[dim_c]] <- fd_temp
        fd_chem[[dim_c]]$fdnames$reps <- chemical_labels[load_number, 2]
      } else {
        fd_chem[[dim_c]]$coefs <- cbind(fd_chem[[dim_c]]$coefs, fd_temp$coefs)
        fd_chem[[dim_c]]$fdnames$reps <-
          append(
            fd_chem[[dim_c]]$fdnames$reps,
            chemical_labels[load_number, 2]
          )
      }
    }
  }

  setwd("../")
  return(list(
    fd = fd_chem, groupd = chemical_labels$Class[keep_ind],
    rep_names = fd_chem[[1]]$fdnames$reps
  ))
}


models <- c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk")

hddc_chem <- function(init = "random", nbasis = 10) {
  trig <- merge_chemical_csv(nbasis = nbasis, allowed_classes = 1:10)
  thresh <- c(0.01, 0.05, 0.1, 0.2)
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd == 1 | trig$groupd == 10,
    as.character(trig$groupd), "other"
  )
  print(paste0(init, " Starting Point"))
  for (t in thresh) {
    print("----------------------------------------------")
    print(paste0("Threshold = ", t))
    res <- funHDDC(trig$fd,
      K = 3, model = models, threshold = t, itermax = 200,
      nb.rep = 50, init = init
    )
    print(table(trig$groupd, res$class))
    print(table(trig$class, res$class))
    print(res$d)
    for (i in unique(trig$groupd)) {
      for (j in unique(res$class)) {
        print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
        print(as.character(trig$rep_names)[trig$groupd == i & res$class == j])
      }
    }
    print("----------------------------------------------")
  }
}

thddc_chem <- function(init = "random", nbasis = 10) {
  trig <- merge_chemical_csv(nbasis = nbasis, allowed_classes = 1:10)
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd == 1 | trig$groupd == 10,
    as.character(trig$groupd), "other"
  )
  thresh <- c(0.01, 0.05, 0.1, 0.2)
  dconstrain <- c("yes", "no")
  print(paste0(init, " Starting Point"))
  for (t in thresh) {
    for (dcon in dconstrain) {
      print("----------------------------------------------")
      print(paste0("Threshold = ", t, " dconstr = ", dcon))
      res <- tfunHDDC(trig$fd,
        K = 3, model = models, threshold = t, itermax = 200,
        nb.rep = 50, init = init, dfstart = 50, dfupdate = "numeric",
        dconstr = dcon
      )
      print(table(trig$groupd, res$class))
      print(table(trig$class, res$class))
      print(res$d)
      print(res$nux)
      for (i in unique(trig$groupd)) {
        for (j in unique(res$class)) {
          print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
          print(as.character(trig$rep_names)[trig$groupd == i & res$class == j])
        }
      }
      print("----------------------------------------------")
    }
  }
}

chddc_chem <- function(init = "random", nbasis = 10) {
  trig <- merge_chemical_csv(nbasis = nbasis, allowed_classes = c(1, 10))
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd == 1 | trig$groupd == 10,
    as.character(trig$groupd), "other"
  )
  thresh <- c(0.01, 0.05, 0.1, 0.2)
  alpha <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  print(paste0(init, " Starting Point"))
  for (t in thresh) {
    for (al in alpha) {
      print("----------------------------------------------")
      print(paste0("Threshold = ", t, " alphamin = ", al))
      res <- funHDDC1(trig$fd,
        K = 3, model = models, threshold = t, itermax = 200,
        nb.rep = 50, init = init, alphamin = al
      )
      print(table(trig$groupd, res$class))
      print(table(trig$groupd, res$outlier))
      print(table(trig$class, res$class))
      print(res$d)
      print(res$alphax)
      print(res$etax)
      for (i in unique(trig$groupd)) {
        for (j in unique(res$class)) {
          print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
          print(as.character(trig$rep_names)[trig$groupd == i & res$class == j])
        }
      }
      print("----------------------------------------------")
    }
  }
}


plot_chemical_curves <- function(fd_chem) {
  CL <- c("1" = "red", "10" = "blue", "other" = "green")
  fd_chem$groupd <- ifelse(fd_chem$groupd == 1 | fd_chem$groupd == 10,
    as.character(fd_chem$groupd), "other"
  )
  par(mfrow = c(3, 4))
  for (i in 1:length(fd_chem$fd)) {
    plot(fd_chem$fd[[i]], col = CL[fd_chem$groupd])
    title(main = paste0("C", i))
  }
}
