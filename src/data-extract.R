library(dplyr)
library(tidyverse)
library(funHDDC)
library(fda)
library(rworldmap)
library(RColorBrewer)
library(mclust)
#library(segclust2d)
setwd("~/cmpt496-mining-covid19") #Working directory
# add a function to plot the fixed time series data and to plot the unfixed data

CL <- brewer.pal(n = 8, name = "Set1")
CLUSTERIDENTIFIER <- "cluster"

fix_negatives <- function(fix_data, column_to_fix) {
  #########################################################
  # Input: fix_data(dataframe) -> dataframe to fix        #
  #        column_to_fix(string) -> column to fix         #
  # Purpose: Assumes that the column to fix is some 'new' #
  #  column with an available 'total' column. Fixes all   #
  #  negative values in the given column.                 #
  # Output: (dataframe) res -> data formatted into form   #
  #                                                       #
  #########################################################
  
  # find coulmns that need to be fixed and assign fix ratios
  tofix <- filter(fix_data, !!as.symbol(column_to_fix)<0)
  tofix$ratio <- 1+dplyr::select(tofix, !!as.symbol(column_to_fix))/(dplyr::select(tofix, !!as.symbol(gsub("new", "total", column_to_fix, fixed = TRUE)))-dplyr::select(tofix, as.symbol(column_to_fix)))
  
  for( i in 1:nrow(tofix)){ #apply fix to all previous days
    # retrieve the column to fix
    colval <- tofix[i,]
    # extract data from the column
    startDate <- colval$date
    isoc <- colval$iso_code
    ra <- colval$ratio
    # find first and last index to fix
    s_ind <- head(which(fix_data$iso_code==isoc), 1)
    e_ind <- which(fix_data$iso_code==isoc & fix_data$date==startDate)
    # check the index
    # fix data
    fix_data[(s_ind:e_ind-1),column_to_fix] <- fix_data[(s_ind:e_ind-1),column_to_fix]*as.double(ra)
    fix_data[e_ind,column_to_fix] <- mean(fix_data[((e_ind-7):(e_ind-1)),column_to_fix])
  } 
  return(fix_data)
}

reformat_data <- function(data_to_format, country_data, columns) {
  #########################################################
  # Input: data_to_format(dataframe) -> dataframe to      #
  #         format into a weekly column form              #
  #        country_data(dataframe) -> countries to use    #
  #        columns(vector(string)) -> columns to use      #
  # Purpose: Format data_to_format into a weekly mean     #
  #  form on the columns given                            #
  # Output: (dataframe) fix_data -> data fixed at column  #
  #                                                       #
  #########################################################
  
  data_to_format[is.na(data_to_format)] <- 0
  # get number of weeks
  # maybe change to lowest date to highest date
  nweeks <- as.integer(difftime(max(data_to_format$date), min(data_to_format$date)))%/%7
  
  # make named columns for data
  Cnames <- c("Unnamed","iso_code","country")
  n_cols <- length(columns)
  for(i in 1:nweeks) { #crates col names for data
    for(j in 1:n_cols) {
      Cnames <- c(Cnames,paste0("W",i,"_", columns[j]))
    }
    
  }
  # make this a general function so we can make it with unfixed data too
  n_rows <- nrow(country_data)
  res <- data.frame(matrix(0,ncol = 3+4*nweeks, nrow = n_rows)) #dataframe for weekly numbers
  colnames(res) <- Cnames
  for (i in 1:n_rows){
    bbb <- country_data[i,2]
    aaa <- data_to_format[data_to_format$iso_code==bbb,]
    res[i,1] <- i
    res[i,2] <- bbb
    res[i,3] <- country_data[i,1]
    
    for(j in 1:nweeks){
      # min(f1$date)
      lll <- (as.Date(min(data_to_format$date))+(j-1)*7)  # THIS HAS BEEN MODIFIED CHECK TO MAKE SURE IT IS OK
      uuu <- (as.Date(min(data_to_format$date))+(j-1)*7+6)
      ccc <- aaa[lll<=aaa$date & aaa$date<=uuu,] #select data for week
      if (dim(ccc)[1]>0){
        for(k in 1:n_cols) {
          nc <- mean(ccc[,c(columns[k])]) #fill in number of cases column
          res[i,4 + n_cols*(j-1)+(k-1)] <- nc
        }
        
      }
    }
  }
  return(res)
  
}

gen_CovidData <- function(filename, raw = FALSE) {
  #########################################################
  # Input: None                                           #
  # Purpose: Extract and reformat data from owid data     #               
  # Output: (dataframe) res -> reformated data            #
  #                                                       #
  #########################################################
  
  # extract owid data
  f1 <- read.csv(file="./owid-data/owid-covid-data.csv") #this is the input file from OWID
  f1 <- f1[c("iso_code",	"continent","location","date","new_cases_per_million","new_deaths_per_million","total_cases_per_million","total_deaths_per_million", "stringency_index", "new_tests_per_thousand")]
  
  # trim unneeded information here, we can use this as our filter with another function
  f2 <- read.csv(file="SELECTED.csv") #this are the selected countries we work with
  
  # merge along country code and iso code
  f1 <- merge(f1,f2, by.y="Country.Code", by.x="iso_code")

  f1[is.na(f1)] <- 0
  
  # fix negatives in the columns that we need to
  if(!raw) {
    f1 <- fix_negatives(fix_data = f1, column_to_fix = "new_cases_per_million")
    f1 <- fix_negatives(fix_data = f1, column_to_fix = "new_deaths_per_million")
  }
  
  # get the reformatted data
  res <- reformat_data(data_to_format = f1, country_data = f2, columns = c("new_cases_per_million", "new_deaths_per_million", "new_tests_per_thousand", "stringency_index"))
  res[is.na(res)] <- 0
  write.csv(res, file = paste0(filename, ".csv"))
  return(res)
}

# play with the options in this
cluster <- function(dtc, nclusters = 2:6, t = 0.2, wtn = 0) {
  models <- c('AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk')
  # redirect the output of funHDDC clusters_nsplines_threshold.txt
  sink(paste0("k", head(nclusters, 1), "-", tail(nclusters, 1), "_s", dtc$basis$nbasis, "_t", t, "_" , dtc$dataname, "_trial", wtn, ".txt")) # give this a proper name
  cluster <- funHDDC(dtc, K=nclusters, model=models, threshold = t, init='random', nb.rep = 100)
  sink()
  # save information from original dtc in cluster
  cluster$rownames <- dtc$rownames
  cluster$nsplines <- dtc$basis$nbasis
  cluster$weeks <- dtc$weeks
  cluster$dname <- dtc$dataname
  cluster$thresh <- t
  cluster$nclusters <- nclusters
  return(cluster)
  
}

format_fd <- function(df, cc, nsplines = -1, 
                      weeks = 1:ncol(dplyr::select(df, 
                                                   contains(cc, 
                                                            ignore.case = TRUE, 
                                                            vars = NULL))),
                      correction=0, max_ratio=1) {
  # save gcv and plot the results
  dft <- dplyr::select(df, contains(cc, ignore.case = TRUE, vars = NULL))
  if(max(weeks)>ncol(dft)){ # if we go over the max weeks set end to the max
    weeks <- weeks[1]:ncol(dft)
  }
  dft <- dft[weeks]
  m <- data.matrix(t(dft)) + correction
  
  if(nsplines < 0){ # if default number of splines then find the best one
    # initialize checking
    best_gcv <- 0
    best_sbd <- NULL
    gcv_results <- c()
    for(i in seq(4, floor((length(weeks)-1)/max_ratio))){ # avoid last week as it overfits
      # set a basis and check result
      basis <- create.fourier.basis(c(0,1), nbasis=i)
      sbd <- smooth.basis(argvals=seq(0,1,length.out = nrow(m)),y=m,fdParobj=basis)
      if(best_gcv == 0 || colSums(as.matrix(sbd$gcv)) < best_gcv){ # save best result
        best_sbd <- sbd
        best_gcv <- colSums(as.matrix(sbd$gcv))
      }
      # save results to a vector
      gcv_results <- append(gcv_results, colSums(as.matrix(sbd$gcv)))
    }
    sbd <- best_sbd
    # output the vector of results
    png("gcv_results_temp.png")
    plot(gcv_results)
    dev.off()
  }else{ # or else use the given basis
    basis <- create.fourier.basis(c(0,1), nbasis=nsplines)
    sbd <- smooth.basis(argvals=seq(0,1,length.out = nrow(m)),y=m,fdParobj=basis)
  }
  
  fd <- sbd$fd
  # save extra information about row and col names and weeks kept
  fd$rownames <- c(df$iso_code)
  fd$colnames <- 1:ncol(m)
  fd$weeks <- paste0(head(weeks, 1), "-", tail(weeks, 1))
  fd$dataname <- cc
  return(fd)
}

fix_functional_data <- function(fdata, data, cc="cases", min_accept=-1) {
  points <- fdata$basis$params
  removed <- c()
  for(i in points) {
    values <- eval.fd(i, fdata)
    for(j in 1:ncol(values)) {
      if(values[,j] < min_accept) {
        removed <- append(removed, fdata$rownames[j])
      }
    }
  }
  removed <- unique(removed)
  keep <- setdiff(fdata$rownames, removed)
  use_data <- data %>% filter(iso_code %in% keep)
  new_fd <- format_fd(use_data, cc, nsplines=fdata$basis$nbasis)
  return(list(fixed_data=use_data, fixed_fd=new_fd, iso_k=keep, iso_r=removed))
  
}

covid_cluster_nc <- function(cdata, wtn = 0) {
  if(wtn>0){
    return(paste0("k", max(cdata$class), "_t", cdata$thresh, "_s", cdata$nsplines, "_m", cdata$model, "_", cdata$dname, "_", cdata$weeks, "_trial", wtn))
  }else{
    return(paste0("k", max(cdata$class), "_t", cdata$thresh, "_s", cdata$nsplines, "_m", cdata$model, "_", cdata$dname, "_", cdata$weeks))
  }
}

cluster_to_csv <- function(cdata, wtn = 0) {
  iso_codes <- cdata$rownames
  cats <- cdata$class
  # format into dataframe by iso code and group
  cdf <- data.frame(iso_codes, cats)
  
  # number of classes_weeks_cluster.csv
  write.csv(cdf,file=paste0(covid_cluster_nc(cdata, wtn = wtn),"_", CLUSTERIDENTIFIER, ".csv")) #write to file
}

plot_cluster_map <- function(cdata, save = FALSE, wtn = 0) {
  
  # extract isocodes, categories, and # of categories
  if(save) {
    png(paste0(covid_cluster_nc(cdata, wtn = wtn),"_map.png"))
  }
  iso_codes <- cdata$rownames
  cats <- cdata$class
  numcats <- max(cats)
  
  # combine into joinCountryData2Map()
  df <- data.frame(iso_codes, cats)
  mapdata <- joinCountryData2Map(df, joinCode = "ISO3", nameJoinColumn = "iso_codes")
  
  # plot on map
  mapCountryData(mapdata, nameColumnToPlot = "cats", numCats = numcats, 
                 catMethod = "categorical", colourPalette = CL[1:numcats], 
                 mapTitle = paste0("Threshold: ", cdata$thresh, " Cluster Map"))
  if(save) {
    dev.off()
  }  
}

plot_fd_by_iso <- function(fdata, iso, ofdata = NULL, save = FALSE) {
  # try to make ofdata dotted lines and categorize the iso by colour
  pos <- which(fdata$rownames %in% iso)
  cols = CL[1:length(iso)]
  plot.fd(fdata[pos[1],], type = 'l', lty = 1, col = cols[1])
  if(length(iso) > 1) {
    for(i in 2:length(pos)){
      lines(fdata[pos[i],], type = 'l', lty = 1, col = cols[i])
    }
  }
  if(!is.null(ofdata)) {
    posun <- which(ofdata$rownames %in% iso)
    # plot the opposing data next to the otyher data
    # lines(ofdata[posun,], col = cols)
    for(i in 1:length(posun)){
      lines(ofdata[posun[i],], type = 'l', lty = 2, col = cols[i])
    }
    # deparse(substitute(var)) gets the original variable's name that was passed in for labels
    legend("topright", legend = c(deparse(substitute(fdata)), deparse(substitute(ofdata))), lty = 1:2)
  }
  legend("topleft", legend = iso, col = cols, lty = 1)
  
}

plot_fd_cluster_intra <- function(cdata, fdata, save = FALSE, wtn = 0) {
  
  # extract number of classes
  u <- max(cdata$class)
  
  # for each class search for only rows which would be in each class and plot based off the count
  for(i in 1:u){ # for each group plot it
    if(save) {
      png(paste0(covid_cluster_nc(cdata, wtn = wtn),"_intra_graph_", i,".png"))
    }
    ind <- which(cdata$class==i)
    plot.fd(fdata[ind,], col = CL[i])
    title(main = paste0("Group ", i, " Comparison"), xlab = "Time", ylab = cdata$dname)
    if(save) {
      dev.off()
    }
  }
  
}
# make this work on the multivariate case
plot_fd_cluster_inter <- function(cdata, fdata, save = FALSE, wtn = 0) {
  # extract number of classes and class numbers
  cc <- cdata$class
  u <- max(cdata$class)
  col_codes <- CL[1:u]
  if(save) {
    png(paste0(covid_cluster_nc(cdata, wtn = wtn),"_inter_graph.png"))
  }
  
  # for each class find and replace the values with colors
  for(i in 1:u){
    cc <- replace(cc, cc == i, col_codes[i])
  }
  # 2021-09-16 Iain: go through each functional data in fdata if its a list
  plot.fd(fdata, col = cc)
  title(main = paste0("Comparison Between ", u, "Groups"), xlab = "Time", ylab = cdata$dname)
  if(save) {
    dev.off()
  }
}


genAriMatrix <- function(outfname, identifier = CLUSTERIDENTIFIER) {
  temp = list.files(pattern = paste0("*", identifier, ".csv"))
  abcd = list()
  
  for(i in 1:length(temp)) {
    cr = read.csv(file=temp[i])
    crm = as.matrix(cr)
    abcd[[i]]=crm
  }
  
  ariM=matrix(0,nrow=length(temp),ncol=length(temp)) #ARI Matrix
  
  for (i in 1:length(temp)){
    
    for (j in 1:length(temp)){
      
      ariM[i,j]=adjustedRandIndex(as.matrix(abcd[[i]]),as.matrix(abcd[[j]]))
      
    }
    
  }#compute ARIMatrix
  
  write.table(ariM, file=paste0(outfname,'.csv'),col.names=F,row.names=F,sep=',' ) #save matrix for further analysis
  
}

segWeek <- function(cdata, cc, seg.var) {
  # get only data from segvar
  dft <- cdata[cdata$iso_code %in% seg.var,]
  # get the order
  seg.order <- dft$iso_code
  # get only data from a cc
  dft <- dplyr::select(dft, contains(cc, ignore.case = TRUE, vars = NULL))
  # find the longest segment of similar values, choose the shortest between the two
  short_int <- 1:ncol(cdata)
  lmin <- 0
  for(i in 1:length(seg.var)) {
    rel_int <- min(which(rel(dft[i,])==0)):max(which(rel(dft[i,])==0))
    if(length(rel_int) < length(short_int)) {
      lmin <- length(rel_int) - length(short_int) + 1
      short_int <- rel_int
    }
  }
  # if this difference is too small then add some days
  if(lmin <= 2) {
    lmin <- lmin + 2
  }
  # trim both and add to a new dataframe
  formdata <- data.frame(t(dft[-short_int]))
  colnames(formdata) <- seg.order
  
  shift_seg <- segmentation(formdata, Kmax=10, seg.var = seg.var, lmin = lmin, scale.variable = TRUE)
  plot(shift_seg)
    

}

bict <- NULL
bic_thresh_track <- function(clust) {
  if(is.null(bict)) {
    bict <<- data.frame(clust$threshold, clust$BIC)
    colnames(bict) <<- c("threshold", "BIC")
  } else {
    temp <- data.frame(clust$threshold, clust$BIC)
    colnames(temp) <- c("threshold", "BIC")
    bict <<- rbind(bict, temp)
  }
}
bic_thresh_close <- function(ncluster, save = FALSE, plot = TRUE) {
  min_thresh <- min(bict$threshold)
  max_thresh <- max(bict$threshold)
  title_string <- paste0("k", ncluster, "_", min_thresh, "_", max_thresh,"_BICvsThresh")
  if(save) {
    write.csv(bict, file = paste0(title_string, "_data.csv"))
    if(plot) {
      png(paste0(title_string, "_graph.png"))
    }
  }
  
  if(plot) {
    plot(x = bict$threshold, y = bict$BIC)
    title(main = title_string, xlab = "Threshold", ylab = "BIC")
  }
  
  if(save && plot) {
    dev.off()
  }
  bict <<- NULL
  
}


main <- function(){
  d <- read.csv(file="fixed_data_new.csv") #this is the input file from OWID
  # d <- gen_CovidData("fixed_data_new")
   
  formatd <- format_fd(d, "new_cases_per_million")
  # track bic vs thresh, maybe put this inside the formatted data
   for(ii in 3:6){
     for(i in 0:13){
       if(i == 0){
         thresh <- 0.001
       } else if(i<=2){
         thresh <- 0.005*i
       } else if(i<=4){
         thresh <- 0.05*(i-2) 
       } else {
         thresh <- 0.1*(i-3) 
       }
       clust <- cluster(formatd, nclusters = ii, t = thresh)
       plot_cluster_map(clust, save = TRUE)
       plot_fd_cluster_inter(clust, formatd, save = TRUE)
       plot_fd_cluster_intra(clust, formatd, save = TRUE)
       cluster_to_csv(clust)
       bic_thresh_track(clust)
     }
     
     bic_thresh_close(ii, save = TRUE, plot = TRUE)
   }
  
  
}

