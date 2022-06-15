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

cfunHDDC_outliers <- function(labels, outliers) {
    # 0 is the outlier
    class1_outliers <- 0
    class2_outliers <- 0
    for (idx in 1:length(labels)) {
        if (labels[idx] == 1 && outliers[idx] == 0) {
            class1_outliers <- class1_outliers + 1
        }
        if (labels[idx] == 2 && outliers[idx] == 0) {
            class2_outliers <- class2_outliers + 1
        }
    }
    cat("Number of outliers in Class 1:", class1_outliers, 
        "/", length(labels[labels == 1]), "\n")
    cat("Number of outliers in Class 2:", class2_outliers, "/", 
        length(labels[labels == 2]), "\n")
}
