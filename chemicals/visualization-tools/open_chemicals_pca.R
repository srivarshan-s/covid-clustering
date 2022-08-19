library(funHDDC)
library(fda)
library(tidyverse)
library(tourr)
library(detourr)

merge_chemical_csv <- function(csv_range=1:65, dim_use=1:11, max_time=NULL, 
                               allowed_classes=c(1,10), nbasis=10) {
  
  setwd("./Raw_CSV")
  
  chemical_labels <- read_csv("chemical_index.csv", show_col_types = F)
  # set up initial empty dimensional data
  dim_data <- list()
  keep_ind <- which(chemical_labels$Class %in% allowed_classes)
  csv_range <- csv_range[keep_ind]
    
  # get the max dimensions
  for(load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name, show_col_types = F)

    if(is.null(max_time) || max(load_df[,"Time"]) < max_time) {
      max_time <- as.numeric(max(load_df$Time))
    }
  }
  # create a basis function based on the highest time can use
  basis <- create.bspline.basis(rangeval = c(0, max_time), nbasis = nbasis)

  fd_chem <- list()
  for(load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name)
    # get only allowed rows
    load_df <- load_df[load_df$Time <= max_time,]
    
    for(dim_c in dim_use) {
      dim_name <- paste0("C", dim_c)
      fd_temp <- smooth.basis(load_df$Time, y=load_df[[dim_name]], 
                              fdParobj=basis)$fd
      if(length(fd_chem) < which(dim_use==dim_c)) {
        fd_chem[[which(dim_use==dim_c)]] <- fd_temp
        fd_chem[[which(dim_use==dim_c)]]$fdnames$reps <- chemical_labels[load_number,2]
      } else {
        fd_chem[[which(dim_use==dim_c)]]$coefs <- cbind(fd_chem[[which(dim_use==dim_c)]]$coefs, fd_temp$coefs)
        fd_chem[[which(dim_use==dim_c)]]$fdnames$reps <- 
          append(fd_chem[[which(dim_use==dim_c)]]$fdnames$reps, 
                 chemical_labels[load_number,2])
      }
      
    }
    
  }
  
  setwd("../")
  return(list(fd=fd_chem, groupd=chemical_labels$Class[keep_ind], 
              rep_names=fd_chem[[1]]$fdnames$reps))
}


models <- c('AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk')

hddc_chem <- function(init="random", nbasis=10) {
  trig <- merge_chemical_csv(nbasis=nbasis, allowed_classes=1:10, dim_use = c(1, 5, 9))
  thresh <- c(0.01,0.05,0.1,0.2)
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd==1 | trig$groupd==10, 
                        as.character(trig$groupd), "other")
  print(paste0(init, " Starting Point"))
  for(t in thresh) {
    print("----------------------------------------------")
    print(paste0("Threshold = ", t))
    res <- funHDDC(trig$fd, K=3, model=models, threshold = t, itermax = 200, 
                   nb.rep = 50, init = init)
    if(!is.null(res$class)) {
      print(table(trig$groupd, res$class))
      print(table(trig$class, res$class))
      print(res$d)
      for(i in unique(trig$groupd)) {
        for(j in unique(res$class)) {
          print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
          print(as.character(trig$rep_names)[trig$groupd==i & res$class==j])
        }
      }
    } else {
      print("Failed")
    }
      print("----------------------------------------------")
    }
}

thddc_chem <- function(init="random", nbasis=10) {
  trig <- merge_chemical_csv(nbasis=nbasis, allowed_classes=1:10, dim_use = c(1, 5, 9))
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd==1 | trig$groupd==10, 
                        as.character(trig$groupd), "other")
  thresh <- c(0.01,0.05,0.1,0.2)
  dconstrain <- c("yes", "no")
  print(paste0(init, " Starting Point"))
  for(t in thresh) {
    for(dcon in dconstrain) {
      print("----------------------------------------------")
      print(paste0("Threshold = ", t, " dconstr = ", dcon))
      res <- tfunHDDC(trig$fd, K=3, model=models, threshold = t, itermax = 200, 
                      nb.rep = 50, init = init, dfstart=50,dfupdate="numeric",
                      dconstr=dcon)
      if(!is.null(res$class)) {
        print(table(trig$groupd, res$class))
        print(table(trig$class, res$class))
        print(res$d)
        print(res$nux)
        for(i in unique(trig$groupd)) {
          for(j in unique(res$class)) {
            print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
            print(as.character(trig$rep_names)[trig$groupd==i & res$class==j])
          }
        }
      } else {
        print("Failure")
      }
      print("----------------------------------------------")
    }
  }
}

chddc_chem <- function(init="random", nbasis=10) {
  trig <- merge_chemical_csv(nbasis=nbasis, allowed_classes=1:10, dim_use = c(1, 5, 9))
  trig$class <- trig$groupd
  trig$groupd <- ifelse(trig$groupd==1 | trig$groupd==10, 
                        as.character(trig$groupd), "other")
  thresh <- c(0.01, 0.05, 0.1, 0.2)
  alpha <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  print(paste0(init, " Starting Point"))
  for(t in thresh) {
    for(al in alpha) {
      print("----------------------------------------------")
      print(paste0("Threshold = ", t, " alphamin = ", al))
      res <- funHDDC1(trig$fd, K=3, model=models, threshold = t, itermax = 200, 
                      nb.rep = 50, init = init, alphamin = al)
      if(!is.null(res$class)) {
        print(table(trig$groupd, res$class))
        print(table(trig$groupd, res$outlier))
        print(table(trig$class, res$class))
        print(res$d)
        print(res$alphax)
        print(res$etax)
        for(i in unique(trig$groupd)) {
          for(j in unique(res$class)) {
            print(paste0("Correct Class = ", i, " and Chosen Class = ", j))
            print(as.character(trig$rep_names)[trig$groupd==i & res$class==j])
          }
        }
      }else{
        print("Failure")
      }
      print("----------------------------------------------")
    }
  }
}


plot_chemical_curves <- function(fd_chem) {
  CL <- c("1"="red", "10"="blue", "other"="green")
  fd_chem$groupd <- ifelse(fd_chem$groupd==1 | fd_chem$groupd==10, 
                           as.character(fd_chem$groupd), "other")
  par(mfrow = c(3,4))
  for(i in 1:length(fd_chem$fd)) {
    plot(fd_chem$fd[[i]], col=CL[fd_chem$groupd])
    title(main=paste0("C",i))
  }
}

get_chemical_data <- function(C, nbasis=10) {
  chem <- merge_chemical_csv(allowed_classes = 1:10, nbasis = nbasis)
  chem_df <- data.frame(matrix(nrow = 63, ncol = nbasis+1))
  chem_names <- c()
  tour_dependence <- c()
  for(j in 1:10) {
    chem_df[,j] <- as.numeric(chem$fd[[C]]$coefs[j,])
    chem_names <- append(chem_names, paste0("C", C,".", j))
    tour_dependence <- append(tour_dependence, c(C))
  }
  names(chem_df) <- append(chem_names, "category")
  chem_df$category <- ifelse(chem$groupd==1 | chem$groupd==10, 
                             as.character(chem$groupd), "other")
  return(chem_df)
}

merge_chemical_dataframe <- function(csv_range=1:65, dim_use=1, max_time=NULL, 
                                     allowed_classes=1:11, nbasis=28, 
                                     pca_dim=10) {
  
  setwd("./Raw_CSV")
  
  chemical_labels <- read_csv("chemical_index.csv", show_col_types = F)
  # set up initial empty dimensional data
  dim_data <- list()
  keep_ind <- which(chemical_labels$Class %in% allowed_classes)
  csv_range <- csv_range[keep_ind]
  
  # get the max dimensions
  for(load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name, show_col_types = F)
    
    if(is.null(max_time) || max(load_df[,"Time"]) < max_time) {
      max_time <- as.numeric(max(load_df$Time))
    }
  }
  # create a basis function based on the highest time can use
  basis <- create.bspline.basis(rangeval = c(0, max_time), nbasis = nbasis)
  eval_data <- matrix(nrow = length(csv_range), ncol = 60)
  
  fd_chem <- list()
  for(load_number in csv_range) {
    load_name <- paste0(load_number, ".csv")
    load_df <- read_csv(load_name, show_col_types = F)
    # get only allowed rows
    load_df <- load_df[load_df$Time <= max_time,]
    
    for(dim_c in dim_use) {
      dim_name <- paste0("C", dim_c)
      fd_temp <- smooth.basis(load_df$Time, y=load_df[[dim_name]], 
                              fdParobj=basis)$fd
      # add evaluate and extract data at certain points
      eval_data[which(csv_range==load_number),] <- as.numeric(eval.fd(seq(0, max_time, 
                                                        length.out=60), 
                                                    fd_temp))
    }
    
  }
  
  setwd("../")
  groupd <- chemical_labels$Class[keep_ind]
  dpca <- prcomp(eval_data, rank.=pca_dim)
  dpca_x <- data.frame(dpca$x)
  dpca_x$category <- ifelse(groupd==1 | groupd==10, 
                            as.character(groupd), "other")
  return(list(dpca_x=dpca_x, data=eval_data, groupd=groupd))
}

# ndim should be 2 or 3
plot_chemical_interactive <- function(chem_df, ndim = 3) {
  detour(chem_df[,1:11], tour_aes(projection=-category, colour=category)) |>
    # tour_path(guided_tour(lda_pp(chem_df$category),d=ndim), fps=60, max_bases = 20) |>
    tour_path(guided_tour(holes()), fps=60, max_bases = 20) |>
    show_scatter(scale_factor = 0.1)
}

# example of how to run
con <- merge_chemical_dataframe(dim_use = 11)
plot_chemical_interactive(con$dpca_x)



# sink()
# print("HDDC")
# sink("hddc_10spline_ldapp_con.txt")
# hddc_chem(init="kmeans", nbasis=10)
# sink()
# print("THDDC")
# sink("thddc_10spline_ldapp_con.txt")
# thddc_chem(init="kmeans", nbasis=10)
# sink()
# print("CHDDC")
# sink("chddc_10spline_ldapp_con.txt")
# chddc_chem(init="kmeans", nbasis=10)
# sink()

