set_seed <- function() {
    set.seed(14071)
}



create_labels <- function(label1, label2) {
    labels <- c()
    for (i in 1:50) {
        labels <- append(labels, label1)
    }
    for (i in 1:40) {
        labels <- append(labels, label2)
    }
    return(labels)
}



functional_data <- function(data, type, nbasis) {
    if (type == "fourier") {
        set_seed()
        basis <- fda::create.fourier.basis(
            rangeval = c(0, 240),
            nbasis = nbasis
        )
        ecg_fdata <- smooth.basis(
            argvals = seq(0, 240, length.out = 241),
            y = t(data.matrix(data)),
            fdParobj = basis
        )$fd
        return(ecg_fdata)
    } else if (type == "bspline") {
        set_seed()
        basis <- create.bspline.basis(
            rangeval = c(0, 480),
            nbasis = nbasis
        )
        ecg_fdata <- smooth.basis(
            argvals = seq(0, 480, length.out = 241),
            y = t(data.matrix(data)),
            fdParobj = basis
        )$fd
        return(ecg_fdata)
    } else {
        stop("Pass valid arguements!")
    }
}



new_functional_data <- function(data, type, nbasis) {
  if (type == "fourier") {
    set_seed()
    basis <- fda::create.fourier.basis(
      rangeval = c(0, 240),
      nbasis = nbasis
    )
    ecg_fdata <- smooth.basis(
      argvals = seq(0, 240, length.out = 241),
      y = t(data.matrix(data)),
      fdParobj = basis
    )$fd
    return(ecg_fdata)
  } else if (type == "bspline") {
    set_seed()
    basis <- create.bspline.basis(
      rangeval = c(0, 240),
      nbasis = nbasis
    )
    ecg_fdata <- smooth.basis(
      argvals = seq(0, 240, length.out = 241),
      y = t(data.matrix(data)),
      fdParobj = basis
    )$fd
    return(ecg_fdata)
  } else {
    stop("Pass valid arguements!")
  }
}



change_labels <- function(labels) {
    for (idx in seq_len(length(labels))) {
        if (labels[idx] == 1) {
            labels[idx] <- 2
        } else {
            labels[idx] <- 1
        }
    }
    return(labels)
}



find_ccr <- function(labels, pred) {
    cf_matrix <- table(labels, pred)
    ccr <- (cf_matrix[1, 1] + cf_matrix[2, 2] + cf_matrix[3, 3])
    ccr <- ccr / sum(cf_matrix)
    return(ccr)
}



find_optimal_labels <- function(labels, pred) {
    max_ccr <- 0
    new_vals <- c(1, 2, 3)
    max_vals <- new_vals
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    new_vals <- c(1, 3, 2)
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    new_vals <- c(2, 1, 3)
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    new_vals <- c(2, 3, 1)
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    new_vals <- c(3, 1, 2)
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    new_vals <- c(3, 2, 1)
    pred_temp <- new_vals[as.factor(pred)]
    ccr <- find_ccr(labels, pred_temp)
    print(ccr)
    if (ccr > max_ccr) {
        max_vals <- new_vals
        max_ccr <- ccr
    }
    return(max_vals[as.factor(pred)])
}



find_misclassified_labels <- function(result, labels, outlier_labels) {
    misclassified_labels <- c()
    for (idx in seq_len(length(labels))) {
        if (labels[idx] != result$class[idx]) {
            misclassified_labels <- append(
                misclassified_labels,
                paste("rep", as.character(idx), sep = "")
            )
        }
    }
    cat("Number of misclassified labels:", length(misclassified_labels), "\n")
    cat(
        "Number of outliers that are misclassified:",
        length(intersect(outlier_labels, misclassified_labels)), "\n"
    )
}

find_eta_values <- function(result) {
    eta_vals <- c(result$etax[1], result$etax[2], result$etax[3])
    cat("Eta values:", eta_vals, "\n")
}

cfunhddc_outliers <- function(labels, outliers) {
    # 0 is the outlier
    class1_outliers <- 0
    class2_outliers <- 0
    class3_outliers <- 0
    for (idx in seq_len(length(labels))) {
        if (labels[idx] == 1 && outliers[idx] == 0) {
            class1_outliers <- class1_outliers + 1
        }
        if (labels[idx] == 2 && outliers[idx] == 0) {
            class2_outliers <- class2_outliers + 1
        }
        if (labels[idx] == 3 && outliers[idx] == 0) {
            class3_outliers <- class3_outliers + 1
        }
    }
    cat(
        "Number of outliers in Class 1:", class1_outliers,
        "/", length(labels[labels == 1]), "\n"
    )
    cat(
        "Number of outliers in Class 2:", class2_outliers, "/",
        length(labels[labels == 2]), "\n"
    )
    cat(
        "Number of outliers in Class 3:", class3_outliers, "/",
        length(labels[labels == 3]), "\n"
    )
}
