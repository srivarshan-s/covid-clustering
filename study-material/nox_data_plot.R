#********************************* RUNNING SIMS *********************************#
# cluster_testdata(1, functional_nox_data(), K=2, 
#                 fun_names=c("funHDDC","cfunHDDC","CNmixt"), nrep=50, 
#                 init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.5))
#plot_sample_data(ianos$plot[2,])
#plot_contaminated_fd(functional_nox_data(),ianos3$plot[2,],show_contam=T,split_class = T, split_contam = F)
#********************************* REQUIRED LIBRARIES *********************************#
library(stats)
library(funHDDC) # clustering library
library(RColorBrewer) # color pallette library
library(ContaminatedMixt)
library(dplyr)
library(tidyverse)
library(fda)
library(rworldmap)
library(RColorBrewer)
library(mclust)
library(fda.usc)
CL <- brewer.pal(n = 8, name = "Set1") # color to use for graphs


split_data <- function(clx, clo, sh_con, sp_cl, sp_con, f_ind=1:length(clx)) {
  # don't worry about this function
  index_split <- list()
  if(sp_cl) {
    # split data into lists based on class
    for(class in unique(clx)) {
      index_split[[class]] <- split_data(clx[clx==class], 
                                         clo[clx==class],
                                         sh_con, F, sp_con,
                                         f_ind=f_ind[clx==class])
      
    }
    return(flatten(index_split))
  } else if(sp_con) {
    for(contam in unique(clo)) {
      index_split[[contam+1]] <- split_data(clx[clo==contam], 
                                            clo[clo==contam],
                                            sh_con, F, F,
                                            f_ind=f_ind[clo==contam])
    }
    return(flatten(index_split))
  } else if(sh_con) {
    clx[clo==0] <- 3 
    index_split[[1]] <- list(ind=f_ind, col=clx)
  } else {
    index_split[[1]] <- list(ind=f_ind, col=clx)
  }
  return(index_split)
}


plot_contaminated_fd <- function(fdata, fdata1,show_contam=T, split_class=F, 
                                 split_contam=F, 
                                 curves=1:length(fdata1$class)) {
  #############################################################################
  # Input:  fdata(list) -> list that contains functional data to plot
  #          this list should be from the simulate_fd() function
  #         show_contam(bool) : default TRUE -> whether to colour
  #          contaminated curves differently in the plot.
  #         split_class(bool) : default FALSE -> whether to plot classes on
  #          separate graphs.
  #         split_contam(bool) : default FALSE -> whether to plot 
  #          contaminated data separately from the uncontaminated data.
  #         curves(vector) : default 1:length(fdata$simdata$clx) -> which
  #          curves to plot of the given data.
  # Purpose: Plot the functional data from the list in the requested
  #  combination of formatting parameters.
  # Output: None -> Displays graphs
  # Usage: 
  # > Display  all the curves together with differently coloured contam
  # plot_contaminated_fd(fdata, show_contam=T, 
  #                      split_class = F, split_contam = F)
  # > Display all contaminated curves separately and in the colour of their
  # > original group separate from uncontaminated curves
  # plot_contaminated_fd(fdata, show_contam=F, 
  #                      split_class = F, split_contam = T)
  # > Display separate groups with their contaminated data on them individual
  # plot_contaminated_fd(fdata, show_contam=T, 
  #                      split_class = T, split_contam = F)
  #############################################################################
  #plot_contaminated_fd(ianos$data,ianos$plot[2,],show_contam=T,split_class = T, split_contam = F)
  # extract class info
  n <- length(fdata1$class)
  # set up plotting with coloration
  clo <- t(fdata1$contam)
  clx <- fdata1$class
  
  # get how to split the plots
  plot_features <- split_data(clx, clo, show_contam, split_class, 
                              split_contam, f_ind=curves)
  print(plot_features)
  # make the plots split up however the data is
  for(i in 1:length(plot_features)) {
    plot(fdata$fd[plot_features[[i]]$ind], col=CL[plot_features[[i]]$col])
  }
  
}
#********************************* IMPORTED DATA *********************************#
# poblenou data set contains functional data
# in $df it contains days of the week and day.festive / groups and contam
data("poblenou")
nox <- poblenou$nox
# univariate data from 2 groups simulated from funHDDC
# last column is the classifications 100 columns of data
data("trigo")

#********************************* GLOBAL VARIABLES =*********************************#
CL <- brewer.pal(n = 8, name = "Set1") # color to use for graphs

#********************************* SIMULATIONS *********************************#

functional_nox_data <- function(ndata = nox, nsplines = 15) {
  # create basis on 24 hours of day
  basis <- create.bspline.basis(rangeval = seq(0,23, length.out=13), nbasis = nsplines)
 # basis <- create.fourier.basis(rangeval = c(0,23), nbasis = nsplines)
  # smooth data on basis with same limits
  nox_fd <- smooth.basis(argvals = seq(0,23, length.out=24), y = t(as.matrix(ndata$data)), 
                         fdParobj = basis)$fd
  # get classifications for the data, include clusters and contamination
 
  
  nox_class <- ifelse(poblenou$df$day.week %in% c(1,2,3,4,5) , 1, 2)
  nox_class <- ifelse(poblenou$df$day.festive %in% c(1), 2, nox_class)
# daca vreau sa le plotez 1 by 1
#  evalbas <- eval.basis(seq(0, 23, length.out = 24), basis)
 # finaldata <- t(nox_fd$coefs)%*%t(evalbas)
 #  plot(nox[101,])
 #points(0:23, finaldata[101,])
  return(list(fd=nox_fd, 
              class=nox_class, 
              contam=as.numeric(poblenou$df$day.festive)-1))
}

functional_trigo_data <- function(tdata = trigo, nsplines = 50) {
  # create basis between 0 and 1 with nsplines given 
  basis <- create.bspline.basis(rangeval = c(0,1), nbasis = nsplines)
  # store the classification and trim it from the matrix
  classification <- tdata[,101]
  trigo_fd <- nox_fd <- smooth.basis(argvals = seq(0,1, length.out=100), 
                                     y = t(tdata[,1:100]), 
                                     fdParobj = basis)$fd
  return(list(fd=trigo_fd, class=classification, contam=rep(0, 100)))
}

cluster_helper <- function(method, fd, init, nb.rep, threshold,etamax, alphamin, ...) {
  x <- list(...)
  print(method)
  models = c( "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
  switch (method,
          funHDDC = return(funHDDC(fd,K=2,nb.rep=nb.rep, 
                                   init=init, threshold=threshold,model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"))),
          cfunHDDC = return(funHDDC1(fd, K=2, nb.rep=nb.rep, 
                                     init=init, threshold=threshold,model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"), alphamin=alphamin)),
          CNmixt = return(CNmixt(t(fd$coefs),G=2, contamination = T, 
                                 initialization = init,alphamin=alphamin, model = models))
  )
}

get_classification <- function(cclust, cmname) {
  switch (cmname,
          funHDDC = return(cclust$class),
          cfunHDDC = return(cclust$class),
          CNmixt = return(getDetection(cclust)$group)
  )
}
get_d <- function(cclust, cmname) {
  switch (cmname,
          funHDDC = return(cclust$d),
          cfunHDDC = return(cclust$d),
          CNmixt = return(0)
  )
}
get_eta <- function(cclust, cmname) {
  switch (cmname,
          funHDDC = return(0),
          cfunHDDC = return(cclust$etax),
          CNmixt = return(getPar(cclust)$eta)
  )
}
get_alpha <- function(cclust, cmname) {
  switch (cmname,
          funHDDC = return(0),
          cfunHDDC = return(cclust$alphax),
          CNmixt = return(getPar(cclust)$alpha)
  )
}
get_outliers <- function(cclust, cmname, ncurves=0) {
  switch (cmname,
          funHDDC = return(rep(NA, ncurves)),
          cfunHDDC = return(cclust$outlier),
          CNmixt = return(ifelse(getDetection(cclust)$innergroup=="*",
                                 1, 0))
          # CNmixt = return(getDetection(cclust)$innergroup)
  )
}

cluster_testdata <- function(ntests, simdata, K=2, clust_FUN=cluster_helper, 
                             fun_names=c("funHDDC","cfunHDDC","CNmixt"), 
                             rseeds=.Random.seed[1:(ntests*2)],
                             nrep=20,init="kmeans", simpre="premade", 
                             threshold=0.2,  alphamin=c(0.5)) {
  ##############################################################################
  # Input: ntests(int) -> a number of tests to run on the given data           #
  #        simdata(list) fd,class -> a list that contains functional data and  #
  #         classifications for that data                                      #
  #        K(vector) : default 2 number of groups to classify using            #
  #        clust_FUN(function) : default cluster_helper -> a function used by  #
  #         the routine to cluster                                             #
  #        fun_names(vector) : default c("funHDDC","cfunHDDC","CNmixt")        #
  #         -> vector of string names for clustering methods                   #
  #        rseeds(vector) : default .Random.seed[1:(ntests*2)]                 #
  #         -> vector of seeds to use in tests                                 #
  #        nrep(int) : default 20 -> number of repetitions to cluster with     #
  #        init(string) : default kmeans -> clustering initialization to use   #
  #        simpre(string) : default premade -> special name for identifying    #
  #        threshold(float) : default 0.2 -> threshold for model evaluation    #
  #        etamax(float) : default 100 -> max eta to use when optimizing       #
  #        alphamin(float) : default 0.5 -> min alpha to use when optimizing   #
  # Purpose: Run several clustering routines on data potentially multiple      #
  #  times and compare results to expected clusters and between routines       #
  # Output: classifications(dataframe) -> data frame of all classifications    #
  ##############################################################################
  # set up containers to hold results
  ncurves <- length(simdata$contam)
  classifications <- data.frame(matrix(ncol = 4 + ncurves, nrow=0))
  pclass=list()
  #names(pclass) <- c("meth", "outliers", "class")
  names(classifications) <- c("n", "seed", "method", "ari", 1:ncurves) # rename columns
  # prefix for premade data
  threshold1=threshold*1000
  alphamin11=alphamin[1]*100
  csv_prefix <- paste0("n_", ntests, "_K_", min(K), "-", max(K), 
                       "_nm_", length(fun_names), "_nrep_", nrep,"_thresh_",
                       threshold1,"_alpha_", alphamin11,"_", simpre, "_")
  for(i in 1:ntests) {
    current_class <- data.frame(matrix(ncol = 4 + ncurves, nrow=0))
    current_pclass=list()
    oclassif <- list()
    cclassif <- list()
    df<- list()
    deta<- list()
    dalpha<- list()
    seed_fds <- abs(rseeds[i])
    set.seed(seed_fds)
    # save the grouping data
    current_class <- rbind(current_class, c(i, seed_fds, "true class", NA,
                                            simdata$class))
    current_pclass<-rbind(current_pclass,list(meth="true class", contam=NA,class=simdata$class))
    # set the seed that will be used for clustering
    seed_clust <- abs(rseeds[i+ntests])
    for(j in 1:length(fun_names)) { # get clustering results for each method
      set.seed(seed_clust) # reset the seed
      # get clustering results
      cclust <- clust_FUN(fun_names[j], fd=simdata$fd, 
                          alphamin=alphamin, K=K, nb.rep=nrep, 
                          init=init, threshold=threshold)
      
      if(!is.null(get_classification(cclust, fun_names[j]))) { # check if the clustering was successful
        # get the individual classifications of the data
        outliers <- get_outliers(cclust, fun_names[j], ncurves)
        df[[j]]<-get_d(cclust, fun_names[j])
        deta[[j]]<-get_eta(cclust, fun_names[j])
        dalpha[[j]]<-get_alpha(cclust, fun_names[j])
        oclassif[[j]] <- outliers
        cclassif[[j]] <- get_classification(cclust, fun_names[j])
        ari <- adjustedRandIndex(cclassif[[j]],simdata$class)
        print(paste0("ari=",ari ))
        print(paste0("d=",df[[j]] ))
        print(paste0("eta=",deta[[j]] ))
        print(paste0("alpha=",dalpha[[j]] ))
        current_class <- rbind(current_class, 
                               c(i, seed_clust, 
                                 fun_names[j], ari, cclassif[[j]]))
        current_pclass <- rbind(current_pclass, 
                                list(meth=fun_names[j], contam=oclassif[[j]],class=cclassif[[j]]))
        # save the results of ARI in a list
        
        table_track1 <- table(cclassif[[j]], simdata$class)
        print(table_track1)
        #table_track2 <- table(oclassif[[j]], simdata$contam)
        table_track2 <- table(oclassif[[j]], simdata$class)
        print(table_track2)
      } else {
        print("--- FAILURE ---")
        oclassif[[j]] <- rep(NA, ncurves)
        cclassif[[j]] <- rep(NA, ncurves)
        # save failures as NA
        current_class <- rbind(current_class, 
                               c(i, seed_clust, fun_names[j], 0, 
                                 rep(NA, ncurves)))
        current_pclass <- rbind(current_pclass, 
                                list(meth=fun_names[j],contam=rep(NA, ncurves), 
                                     class=rep(NA, ncurves)))
      } # end of if cclust null
    } # end of for each clust_FUN
    # compare all the outliers data classifications in a table
    oclass_table <- matrix(nrow = length(fun_names), ncol = length(fun_names),
                           dimnames = list(fun_names, fun_names))
    cclass_table <- matrix(nrow = length(fun_names), ncol = length(fun_names),
                           dimnames = list(fun_names, fun_names))
    for(i in 1:length(fun_names))
      for(i in 1:length(fun_names)) {
        oclass_table[i, j] <- adjustedRandIndex(oclassif[[i]], oclassif[[j]])
        cclass_table[i, j] <- adjustedRandIndex(cclassif[[i]], cclassif[[j]])
      }
    write.csv(oclass_table, file = paste0(csv_prefix, "outlier_comparison.csv"))
    write.csv(cclass_table, file = paste0(csv_prefix, 
                                          "classification_comparison.csv"))
    # save the current class as its own file
    names(current_class) <- names(classifications) # rename columns
    write.csv(current_class, file = paste0(csv_prefix, "_test_", i, "_clusters.csv"))
    # add the current cluster to tracker
    classifications <- rbind(classifications, current_class)
    pclass<-rbind(pclass, current_pclass)
  } # end of for each ntest
  # save the results
  write.csv(classifications, file = paste0(csv_prefix, "clusters.csv"))
  # return the results
  res=list(class=classifications, plot=pclass)
  return(res)
}

#********************************* PLOTTING FUNCTIONS *********************************#
plot_sample_data <- function(fund=NULL){
  # get the data from poblenou
  if(is.null(fund)) {
    fund <- functional_nox_data()
  }
  # get contaminated groupings
  d=functional_nox_data()
  func <- ifelse(fund$contam == 1, fund$class, max(fund$class)+1)
  func1 <- ifelse(fund$contam == 1, fund$class, fund$class)
  print(d$fd$fdnames$reps[which(fund$contam==0)])
  # plot every d=functional_nox_data()
  plot.fd(d$fd, col=CL[func1], ylab=fund$meth)
  plot.fd(d$fd, col=CL[func], ylab=fund$meth)
  # plot individual groups
  for(cl in unique(func)){
    plot.fd(d$fd[func==cl], col=CL[cl],ylab=fund$meth)}
  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow=c(2,1))
  for(cl in unique(func1)){
    plot.fd(d$fd[func1==cl], col=CL[cl],ylab=fund$meth)
    
    }
 # col1<-ifelse((fund$class-fund$contam) == 0, fund$class, max(fund$class)+1)
  #col2<-ifelse((2*fund$class-fund$contam) == 3, fund$class, max(fund$class)+1)
  #plot.fd(d$fd[func1==1], col=CL[col1],ylab=fund$meth)
  #plot.fd(d$fd[func1==2], col=CL[col2],ylab=fund$meth)
  
}
funHDDC1  <-
  function(data, K=1:10, model="AkjBkQkDk", alphamin=c(0.5),threshold=0.1, itermax=200, eps=1e-6,init='random', criterion="bic",
           algo='EM', d_select="Cattell",  init.vector=NULL, show=TRUE, mini.nb=c(5, 10),
           min.individuals=2, mc.cores=1, nb.rep=2, keepAllRes=TRUE, kmeans.control = list(), d_max=100){
    
    #Options removed from call
    com_dim=NULL
    scaling=FALSE
    noise.ctrl=1e-8
    
    #nb.rep parameter
    if ((init=="random")&(nb.rep<20)){
      nb.rep=20
    }
    #
    # CONTROLS
    #
    
    call = match.call()  #A@5 macthes values to prameters
    
    .hddc_control(call)
    # Control of match.args:
    criterion = .myAlerts(criterion, "criterion", "singleCharacterMatch.arg", "HDDC: ", c("bic", "icl"))
    algo = .myAlerts(algo, "algo", "singleCharacterMatch.arg", "HDDC: ", c('EM', 'CEM', 'SEM'))
    d_select = .myAlerts(d_select, "d_select", "singleCharacterMatch.arg", "HDDC: ", c("cattell", "bic"))
    init = .myAlerts(init, "init", "singleCharacterMatch.arg", "HDDC: ", c('random', 'kmeans', 'mini-em', 'param', "vector"))
    # We get the model names, properly ordered
    model = .hdc_getTheModel(model, all2models = TRUE)
    # kmeans controls
    kmeans.control = .default_kmeans_control(kmeans.control)
    
    
    if (scaling) {
      stop('Scaling is not implemented!')
    } else scaling <- NULL
    
    BIC <- ICL <- c()
    fdobj = data  #A@5 we work on univariate  model so for us it is NOT a list - CHECK ALWAYS ELSE
    if (class(fdobj)!='list') {x = t(fdobj$coefs)}
    else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
    p <- ncol(x)
    
    #
    # Preparing the parallel
    #
    
    if(d_select=="bic"){
      # If the dimension selection is done with BIC, we don't care of the threshold
      threshold = "bic"
    }
    
    if(max(table(K))>1) warning("The number of clusters, K, is made unique (repeated values are not tolerated).")
    K = sort(unique(K))
    if(any(K==1)){
      # K=1 => only one model required
      K = K[K!=1]
      addrows = data.frame(model="AKJBKQKDK", K=1, threshold)
    } else {
      addrows = c()
    }
    
    mkt_Expand = expand.grid(model=model, K=K, threshold=threshold)
    mkt_Expand = do.call(rbind, replicate(nb.rep, mkt_Expand, simplify=FALSE))
    mkt_Expand = rbind(addrows, mkt_Expand) #no need for several runs for K==1
    
    model = as.character(mkt_Expand$model)
    K = mkt_Expand$K
    threshold = mkt_Expand$threshold
    
    # We transform it into an univariate form
    mkt_univariate = apply(mkt_Expand, 1, paste, collapse= "_")
    
    # Mon 'caller' for LTBM that will be used in multi-cores lapply
    hddcWrapper = function(mkt_univariate, ...){
      
      mkt_splitted =  strsplit(mkt_univariate, "_")
      
      # on retrouve model, K and threshold
      model = sapply(mkt_splitted, function(x) x[1])
      K = sapply(mkt_splitted, function(x) as.numeric(x[2]))
      threshold = sapply(mkt_splitted, function(x) ifelse(x[3]=="bic","bic",as.numeric(x[3])))
      
      # (::: is needed for windows multicore)
      res = "unknown error"
      # try(res <- HDclassif:::funhddc_main(model=model, K=K, threshold=threshold, ...))
      try(res <- .funhddc_main1(model=model, K=K, alphamin=alphamin, threshold=threshold, itermax=itermax, ...),silent=TRUE)
      res
    }
    
    # We reset the number of cores to use
    nRuns = length(mkt_univariate)
    if(nRuns<mc.cores) mc.cores = nRuns
    
    # We swicth to the right number of cores + a warning if necessary
    max_nb_of_cores = parallel::detectCores()
    if(mc.cores>max_nb_of_cores){
      warning("The argument mc.cores is greater than its maximun.\nmc.cores was set to ", max_nb_of_cores)
      mc.cores = max_nb_of_cores
    }
    
    
    #
    # Parallel estimations
    #
    
    if(mc.cores == 1){
      # If there is no need for parallel, we just use lapply // in order to have the same output
      
      par.output = lapply(mkt_univariate, hddcWrapper, fdobj=fdobj, method=d_select, algo=algo, eps=eps, init=init, init.vector=init.vector, mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max)
      #A@5  here mkt_univariate is applied to to fdobj
      
    } else if(Sys.info()[['sysname']] == 'Windows'){
      # we use parLapply:
      
      ## create clusters
      cl = parallel::makeCluster(mc.cores)
      ## load the packages
      loadMyPackages = function(x){
        # loadMyPackages = function(none, myFuns){
        # we load the package
        library(funHDDC)
        # we add the functions in the global env
        # for(i in 1:length(myFuns)){
        # 	funName = names(myFuns)[i]
        # 	fun = myFuns[[i]]
        # 	assign(funName, fun, .GlobalEnv)
        # }
      }
      ## Create the functions to export
      # myFuns = list(hddc_main = hddc_main, hddc_e_step=hddc_e_step, hddc_m_step=hddc_m_step, hdc_getComplexity=hdc_getComplexity, hdc_myEigen=hdc_myEigen, hdclassif_dim_choice=hdclassif_dim_choice, hdclassif_bic=hdclassif_bic)
      # par.setup = parallel::parLapply(cl, 1:length(cl), loadMyPackages, myFuns=myFuns)
      par.setup = parallel::parLapply(cl, 1:length(cl), loadMyPackages)
      
      ## run the parallel
      par.output = NULL
      try(par.output <- parallel::parLapply(cl, mkt_univariate, hddcWrapper, DATA=fdobj, method=d_select, algo=algo, itermax=itermax, eps=eps, init=init, init.vector=ifelse(missing(init.vector), NA, init.vector), mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max))
      
      ## Stop the clusters
      parallel::stopCluster(cl)
      
      if(is.null(par.output)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
      
    } else {
      # we use mclapply
      
      par.output = NULL
      try(par.output <- parallel::mclapply(mkt_univariate, hddcWrapper, DATA=fdobj, method=d_select, algo=algo, itermax=itermax, eps=eps, init=init, init.vector=init.vector, mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max, mc.cores=mc.cores))
      
      if(is.null(par.output)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
      
    }
    
    #
    # The results are retrieved
    #
    
    getElement = function(x, what, valueIfNull = -Inf){
      # attention si x est le modele nul
      if(length(x)==1) return(valueIfNull)
      if(!is.list(x) && !what %in% names(x)) return(NA)
      x[[what]][length(x[[what]])]
    }
    
    getComment = function(x){
      # we get the error message
      if(length(x)==1) return(x)
      return("")
    }
    
    # All likelihoods
    LL_all = sapply(par.output, getElement, what="loglik")
    comment_all = sapply(par.output, getComment)
    
    # If no model is valid => problem
    if(all(!is.finite(LL_all))){
      warning("All models diverged.")
      allCriteria = data.frame(model=model, K=K, threshold=threshold, LL = LL_all, BIC=NA, comment=comment_all)
      res = list()
      res$allCriteria = allCriteria
      return(res)
    }
    
    # We select, for each (Q,K), the best run
    n = nrow(mkt_Expand)
    modelKeep = sapply(unique(mkt_univariate), function(x) (1:n)[mkt_univariate==x][which.max(LL_all[mkt_univariate==x])])
    # => we select only the best models
    LL_all = LL_all[modelKeep]
    comment_all = comment_all[modelKeep]
    par.output = par.output[modelKeep]
    BIC = sapply(par.output, getElement, what="BIC")
    ICL = sapply(par.output, getElement, what="ICL")
    comp_all = sapply(par.output, getElement, what="complexity", valueIfNull=NA)
    model = model[modelKeep]
    threshold = threshold[modelKeep]
    K = K[modelKeep]
    
    # We define the criterion of model selection
    CRIT = switch(criterion,
                  bic = BIC,
                  icl = ICL)
    
    # The order of the results
    myOrder = order(CRIT, decreasing = TRUE)
    
    # On sauvegarde le bon modele + creation de l'output
    qui = which.max(CRIT)
    
    prms = par.output[[qui]]
    prms$criterion = CRIT[qui]
    names(prms$criterion) = criterion
    
    # Other output
    prms$call = call
    # We add the complexity
    names(comp_all) = mkt_univariate[modelKeep]
    prms$complexity_allModels = comp_all
    
    # Display
    if(show){
      if(n>1) cat("cfunHDDC: \n")
      
      model2print = sapply(model, function(x) sprintf("%*s", max(nchar(model)), x))
      K2print = as.character(K)
      K2print = sapply(K2print, function(x) sprintf("%*s", max(nchar(K2print)), x))
      thresh2print = as.character(threshold)
      thresh_width = max(nchar(thresh2print))
      thresh2print = sapply(thresh2print, function(x) sprintf("%s%s", x, paste0(rep("0", thresh_width - nchar(x)), collapse="") ))
      
      # on cree une data.frame
      myResMat = cbind(model2print[myOrder], K2print[myOrder], thresh2print[myOrder],.addCommas(prms$complexity_allModels[myOrder]), .addCommas(CRIT[myOrder]), comment_all[myOrder])
      
      myResMat = as.data.frame(myResMat)
      names(myResMat) = c("model", "K", "threshold", "complexity",toupper(criterion), "comment")
      row.names(myResMat) = 1:nrow(myResMat)
      
      # if no problem => no comment
      if(all(comment_all == "")) myResMat$comment = NULL
      
      print(myResMat)
      
      msg = switch(criterion, bic="BIC", icl="ICL")
      cat("\nSELECTED: model ", prms$model, " with ", prms$K, " clusters.\n")
      cat("Selection Criterion: ", msg, ".\n", sep="")
      
    }
    
    # We also add the matrix of all criteria
    allCriteria = data.frame(model=model[myOrder], K=K[myOrder], threshold=threshold[myOrder], LL=LL_all[myOrder],complexity=prms$complexity_allModels[myOrder], BIC=BIC[myOrder], ICL=ICL[myOrder], rank = 1:length(myOrder))
    
    # we add the comments if necessary
    if(any(comment_all != "")) allCriteria$comment = comment_all[myOrder]
    prms$allCriteria = allCriteria
    
    # If all the results are kept
    if(keepAllRes){
      all_results = par.output
      names(all_results) = mkt_univariate[modelKeep]
      prms$all_results = all_results
    }
    
    # Other stuff
    prms$scaling <- scaling
    prms$threshold <- threshold[qui]
    
    return(prms)
  }

.funhddc_main1 <- function(fdobj, K, model,alphamin, itermax=200,threshold, method, algo, eps, init, init.vector, mini.nb, min.individuals, noise.ctrl, com_dim=NULL, kmeans.control, d_max, ...){
  l1=length(alphamin)
  if (l1<K) {alphamin<-c(alphamin,rep(alphamin[1],K-l1))}
  ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD")
  if (class(fdobj)!='list') {DATA = t(fdobj$coefs)}
  else {DATA = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) DATA = cbind(DATA,t(fdobj[[i]]$coefs))}
  p <- ncol(DATA)
  N <- nrow(DATA)
  com_ev <- NULL
  
  # We set d_max to a proper value
  d_max = min(N, p, d_max)
  
  
  if (K>1){
    t <- matrix(0, N, K)
    
    ##############################################
    #initialization for t  
    if(init == "vector"){
      init.vector = unclass(init.vector)
      name <- unique(init.vector)
      for (i in 1:K) t[which(init.vector==name[i]), i] <- 1
    } else if (init=='kmeans') {
      kmc = kmeans.control
      cluster <- kmeans(DATA, K, iter.max=kmc$iter.max, nstart=kmc$nstart, algorithm=kmc$algorithm, trace=kmc$trace)$cluster
      for (i in 1:K) 
        t[which(cluster==i), i] <- 1
      #A@5 initialize t with the result of k-means
    } else if (init=='mini-em'){
      prms_best <- 1
      for (i in 1:mini.nb[1]){
        prms <- .funhddc_main1(DATA, K, model,alphamin, etamax, threshold, method, algo, mini.nb[2], 0, 'random', mini.nb = mini.nb, min.individuals = min.individuals, noise.ctrl = noise.ctrl, com_dim = com_dim, d_max=d_max)
        if(length(prms)!=1){
          if (length(prms_best)==1) prms_best <- prms
          else if (prms_best$loglik[length(prms_best$loglik)]<prms$loglik[length(prms$loglik)]) prms_best <- prms
        }
      }
      
      if (length(prms_best)==1) return(1)
      t <- prms_best$posterior
    } else { #A@5 INIT IS RANDOM  
      t <- t(rmultinom(N, 1, rep(1/K, K))) #A@5 some multinomial
      compteur=1
      while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N, 1, rep(1/K, K)))
      if(min(colSums(t))<1) return("Random initialization failed (n too small)")
    }
  } else t <- matrix(1, N, 1) #@5 K IS 1
  
  # 
  # A@5 initialize v(i)=0.99 - it's a matrix (N by K)
  # A@5 initialize eta(k)=1.01 - it's a vector (K by 1)
  # A@5 
  #A@5 initialization fo alpha
  vx <- matrix(0.99, N, K)
  alphax <- rep(0.999,K)
  etax   <- rep(1.01,K)
  #etax   <- rep(1,K)
  #A@5 initialization for v
  #vx[which(cluster==i), i] <- 0.999
  #for (i in 1:K) 
  #vx[which(t[,i]>0), i] <- 0.999
  #vX   <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  
  
  
  likely <- c()
  I <- 0
  test <- Inf
  while ((I <- I+1)<=itermax && test>=eps){
    #  while ((I <- I+1)<=itermax ){ 
    # print(paste0("Iter=", I))
    
    # A@5 loops is here until itermax or test <eps
    
    if (algo!='EM' && I!=1) t <- t2
    
    # A@5 t2 is probably the new t which gets transfered into t for CEM and SEM
    # Error catching
    if (K>1){
      if(any(is.na(t))) return("unknown error: NA in t_ik")
      
      if(any(colSums(t>1/K)<min.individuals)) return("pop<min.individuals")
    }
    
    m <- .funhddc_m_step1(fdobj, K, t,vx,alphax,etax, model,alphamin,threshold, method, noise.ctrl, com_dim, d_max) # A@5mstep is applied here
    #@A5 here we need to add vi and eta after t
    # print(m$mu)
    
    to <- .funhddc_e_step1(fdobj, m,t) #A@5 estep is applied here
    L <- to$L # A@5 this s the likelyhood  THIS HAS TO BE CHANGED
    t <- to$t#A@5 this is the new t
    vx<-to$vx
    #if (I>3)
    etax<-m$etax
    # for (i in 1:K) if(m$etax[i]<etamax)  etax[i]<-m$etax[i]
    #print(paste0("I=",I, "etax=",etax ))
    #A@5 add v there and return it something like v=t$v;
    
    if (algo=='CEM') { #@A5 FIND ABOUT CEM
      t2 <- matrix(0, N, K)
      t2[cbind(1:N, max.col(t))] <- 1
    } else if(algo=='SEM') {
      t2 <- matrix(0, N, K)
      for (i in 1:N)	t2[i, ] <- t(rmultinom(1, 1, t[i, ]))
    }
    #likelihood contrib=L
    likely[I] <- L
    #if (I!=1) test <- abs(likely[I]-likely[I-1]) # A@5 did it converge? 
    if (I==2) test <- abs(likely[I]-likely[I-1])
    else if (I>2)
    {
      lal<-(likely[I]-likely[I-1])/(likely[I-1]-likely[I-2])
      lbl<-likely[I-1]+(likely[I]-likely[I-1])/(1.0-lal)
      test<-abs(lbl-likely[I-1])
    }
    #  print(paste0("lik=",likely[I] ))
    # print(paste0("alphax=",m$alphax ))
    #  print(paste0("etax=",m$etax ))
  }
  # print(paste0("I=",I,"lik=",likely[I-1] ))
  #  print(paste0("alphax=",m$alphax ))
  # print(paste0("etax=",m$etax ))
  # We retrieve the parameters
  
  #A@5 CAN'T TOUCH THIS
  #--------------------------------------
  # a
  if ( model%in%c('AKBKQKDK', 'AKBQKDK', 'AKBKQKD', 'AKBQKD') ) {
    a <- matrix(m$a[, 1], 1, m$K, dimnames=list(c("Ak:"), 1:m$K))
  } else if(model=='AJBQD') {
    a <- matrix(m$a[1, ], 1, m$d[1], dimnames=list(c('Aj:'), paste('a', 1:m$d[1], sep='')))
  } else if ( model%in%c('ABKQKDK', 'ABQKDK', 'ABKQKD', 'ABQKD', "ABQD") ) {
    a <- matrix(m$a[1], dimnames=list(c('A:'), c('')))
  } else a <- matrix(m$a, m$K, max(m$d), dimnames=list('Class'=1:m$K, paste('a', 1:max(m$d), sep='')))
  
  # b
  if ( model%in%c('AKJBQKDK', 'AKBQKDK', 'ABQKDK', 'AKJBQKD', 'AKBQKD', 'ABQKD', 'AJBQD', "ABQD") ) {
    b <- matrix(m$b[1], dimnames=list(c('B:'), c('')))
  } else b <- matrix(m$b, 1, m$K, dimnames=list(c("Bk:"), 1:m$K))
  
  # d, mu, prop
  d <- matrix(m$d, 1, m$K, dimnames=list(c('dim:'), "Intrinsic dimensions of the classes:"=1:m$K))
  mu <- matrix(m$mu, m$K, p, dimnames=list('Class'=1:m$K, 'Posterior group means:'=paste('V', 1:p, sep='')))
  prop <- matrix(m$prop, 1, m$K, dimnames=list(c(''), 'Posterior probabilities of groups'=1:m$K))
  alphax <- matrix(m$alphax, 1, m$K, dimnames=list(c(''), 'Probabilities of outliers of groups'=1:m$K))
  etax <- matrix(m$etax, 1, m$K, dimnames=list(c(''), 'Inflation coefficient of groups'=1:m$K))
  
  
  # A25 add alpha and eta
  
  # Other elements
  complexity <- .hdc_getComplexity(m, p)
  class(etax)<-class(alphax)<-class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- 'hd' #A@5 alpha and eta will nedd class too
  cls <- max.col(t)#@ here the cluster is found
  outl=matrix(0,N,1)
  for(i in 1:N)
  {
    outl[i]<-(to$vx[i,cls[i]]>0.5)+0
  }
  ## we have to add the clasification as outlier or not
  converged = test<eps
  
  params = list(model=model, K=K, d=d, a=a, b=b, mu=mu, prop=prop,alphax=alphax,etax=etax, ev=m$ev, Q=m$Q,fpca=m$fpcaobj, loglik=likely[length(likely)], loglik_all = likely, posterior=t, class=cls,outlier=outl, com_ev=com_ev, N=N, complexity=complexity, threshold=threshold, d_select=method, converged=converged) #A25 add alpha and eta
  #  params = list(model=model, K=K, d=d, a=a, b=b, mu=mu, prop=prop,alphax=alphax,etax=etax, ev=m$ev, Q=m$Q,fpca=m$fpcaobj, loglik=likely[length(likely)], loglik_all = likely, posterior=t, class=cls, com_ev=com_ev, N=N, complexity=complexity, threshold=threshold, d_select=method, converged=converged) #A25 add alpha and eta
  
  #A@5 CONTINUE FROM HERE
  
  
  # We compute the BIC / ICL
  bic_icl = .hdclassif_bic(params, p)
  params$BIC = bic_icl$bic
  params$ICL = bic_icl$icl
  
  # We set the class
  class(params) <- 'cfunHDDC'
  
  return(params)
}

.funhddc_e_step1  <- function(fdobj, par,t){ #A@5 this is the E step
  if (class(fdobj)!='list') {x = t(fdobj$coefs)}  #A@5 THIS IS OUR CASE
  else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))} #NOT THIS
  p <- ncol(x)
  N <- nrow(x)
  K <- par$K
  a <- par$a
  b <- par$b
  mu <- par$mu
  d <- par$d
  prop <- par$prop
  Q <- par$Q
  alphax<-par$alphax
  etax<-par$etax
  #A@5 add the alhpas and the etas
  
  b[b<1e-6] <- 1e-6
  
  #############################aici am ajuns
  #################################################
  K_pen <- matrix(0,K,N)
  # L1_pen <- matrix(0,K,N)
  #  L2_pen <- matrix(0,K,N)
  S_pen<-matrix(0,K,N)
  Sg_pen<-matrix(0,K,N)
  Sb_pen<-matrix(0,K,N)
  mah_pen<-matrix(0,K,N)
  ax_pen<-matrix(0,N,K)
  #pdf_pen<-matrix(0,N,K)
  Stot_pen<-matrix(0,K,N)
  vx <- matrix(0.99,N,K)
  s=rep(0,K)
  #A@5 here is the work
  for (i in 1:K) { #A@5 this must become formula 12,13,14,15 lot's of work
    s[i] <- sum(log(a[i,1:d[i]])) #A@5 this is the third term in formula 15
    #mah_pen[i,]=imahalanobis(x, i,mu,par$fpcaobj[[i]]$W ,Q, a, d, b)
    Qk=Q[[i]]
    aki=sqrt(diag(1/a[i,1:d[i]],d[i]))
    bki=b[i]
    muki=mu[i,]
    Wki=par$fpcaobj[[i]]$W
    mah_pen[i,]=.imahalanobis(x, muki,Wki ,Qk, aki, bki)
    S_pen[i,]=1/2*(p*log(etax[i])+(1/etax[i]-1)*mah_pen[i,])
    ax_pen[,i]=t(alphax[i]*exp((S_pen[i,]))+(1-alphax[i]))
    
    vx[,i]=t(alphax[i]*exp(S_pen[i,]))/ax_pen[,i]
    K_pen[i,] <- 1.0/etax[i]*mah_pen[i,]+p*log(etax[i])+s[i]+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi) #A@5 this is fomula 15. but without p*log(2*pi)
    ##vx like in CNmixt
    # Sg_pen[i,]=-0.5*(mah_pen[i,]+s[i]+(p-d[i])*log(b[i])+p*log(2*pi))
    #Sb_pen[i,]=-0.5*(1.0/etax[i]*mah_pen[i,]+p*log(etax[i])+s[i]+(p-d[i])*log(b[i])+p*log(2*pi))
    #Stot_pen[i,]=alphax[i]*exp(Sg_pen[i,])+(1-alphax[i])*exp(Sb_pen[i,])
    #vx[,i]=t(alphax[i]*exp(Sg_pen[i,])/Stot_pen[i,])    
    # j11<-quantile(vx[,i], c(1-alphax[i]))
    #vx[vx[,i]>j11,i] <- 0.999
    #vx[which(i!=max.col(t)[]),i]=0.999
  }
  
  
  vx[vx<1e-6] <- 1e-6
  
  #pdf_pen=t(prop*Stot_pen)
  #############################################
  
  
  
  # vx[which(Stot_pen<1e-10)]=0
  
  #what is this L?
  A <- -1/2*t(K_pen) 
  L<-sum(log(rowSums(ax_pen*exp(A-apply(A,1,max))))+apply(A,1,max))
  #likelihood as in CNmixt
  # L<-sum(log(rowSums(pdf_pen)))
  #A25 add calculations for v (formula 13)exp
  
  t <- matrix(0,N,K)
  # for (i in 1:K) t[,i] <- ax_pen[,i]/rowSums(t(exp(-0.5*K_pen+0.5*matrix(K_pen[i,],K,N,byrow = TRUE)))*ax_pen) #A@5 this  1 over is the exponential part on the denominator on formula 12  (the inverse)
  for (i in 1:K) t[,i] <- ax_pen[,i]*exp(A[,i]-apply(A,1,max))/rowSums(ax_pen*exp(A-apply(A,1,max)))
  
  # t like in CN mixt
  # for (i in 1:K) t[,i] <-pdf_pen[,i]/rowSums(pdf_pen)
  #t[t<1e-6] <- 1e-6 
  trow<-numeric(N)
  tcol=numeric(K)
  trow<-rowSums(t)
  tcol=colSums(t)
  if(any(tcol<p)) { #print(paste0("t-problems ",tcol ))
    #for (i in 1:N)
    # t[i,]=(t[i,]+0.0000001)/(trow[i]+K*0.0000001) 
    t=(t+0.0000001)/(trow+K*0.0000001) 
    #print(paste0("t-after ",colSums(t) ))
  }
  list(t=t,vx=vx,L=L) #A@5 add v to t and return it  as well
}

# A@5 pass the z, v, and eta in a names list from the e step par

.funhddc_m_step1  <- function(fdobj, K, t,vx,alphax,etax, model,alphamin, threshold, method, noise.ctrl, com_dim, d_max){ #A@5 this is the M step
  if (class(fdobj)!='list') { x = t(fdobj$coefs) #A@5 THIS IS OUR CASE
  } else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
  #A@5 x is the coefficient in the fdobject
  # Some parameters
  N <- nrow(x)
  p <- ncol(x)
  prop <- c()
  n <- colSums(t) 
  prop <- n/N  #A5 props is pi as a vector
  alphax <- colSums(t*vx)/n
  #print(paste0("calcalphax=",alphax ))
  alphmax <- function(alpha,z, v) sum(z*(v*log(alpha) + (1 - v)*log(1-alpha)))
  for (i in 1:K) {
    if(alphax[i] < alphamin[i]) {
      # z is the vector z[,k] v is a vector v[,k] maximize alpha
      alphax[i] <- optimize(alphmax, c(alphamin[i],1), maximum = TRUE,z = t[,i],v = vx[,i])$maximum
      #  alphax[i]=alphamin[i]
    }
  }
  
  
  correctionX <- array(1,c(N,K),dimnames=list(1:N,paste("comp.",1:K,sep=""))) 
  
  correctionX <- vx+(1.0-vx)*matrix(1.0/etax,N,K,byrow=TRUE)#like fact in CNMixt  
  mu <- matrix(NA, K, p)
  corX=matrix(0,N,K)
  corX=t*correctionX
  for (i in 1:K) {
    #mu[i, ] <- colSums(x*t[, i])/n[i] #A@5 add the big parenthesis from formula 18 and becomes something like colSums(x*t[,i]* a slive of v)
    #A@5 and n[i] must be multiplied by some form of colsum of the product between t and the same thing as above
    # a@5 new
    
    mu[i,]<-apply(t(as.matrix(corX[,i]) %*% matrix(1,1,p)) * t(x), 1, sum) / sum(corX[,i])
    
    #mu[i,]<-colSums(x*matrix(corX[,i],N,p,byrow=FALSE))/sum(corX[,i])
  }
  # print(mu)
  # for (i in 1:K) mu[i, ] <- colSums((z[, i]*(v[, i] + (1-v[, i])/eta[i]))*x)/colSums((z[, i]*(v[, i] + (1-v[, i])/eta[i])))
  ind <- apply(t>0, 2, which)
  n_bis <- c()
  for(i in 1:K) n_bis[i] <- length(ind[[i]])
  
  #
  #Calculation on Var/Covar matrices
  #
  
  # we keep track of the trace (== sum of eigenvalues) to compute the b
  traceVect = c()
  
  
  ev <- matrix(0, K, p)
  Q <- vector(mode='list', length=K) #A@5 this is matrix q_k
  fpcaobj = list()
  for (i in 1:K){
    donnees <- .mypca.fd1(fdobj,t[,i],corX[,i]) #A@5 we must change .mypca- bu everything else stays the same
    #A@% DON"T TOUCH
    #-----------------------------------
    traceVect[i] = sum(diag(donnees$valeurs_propres))
    ev[i, ] <- donnees$valeurs_propres
    Q[[i]] <- donnees$U
    fpcaobj[[i]] = donnees
  }
  
  
  #Intrinsic dimensions selection
  # browser()
  if (model%in%c("AJBQD", "ABQD")){
    d <- rep(com_dim, length=K)
  } else if ( model%in%c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD") ){
    dmax <- min(apply((ev>noise.ctrl)*rep(1:ncol(ev), each=K), 1, which.max))-1
    if(com_dim>dmax) com_dim <- max(dmax, 1)
    d <- rep(com_dim, length=K)
  } else {
    d <- .hdclassif_dim_choice(ev, n, method, threshold, FALSE, noise.ctrl)
  }
  
  #Setup of the Qi matrices
  
  for(i in 1:K) Q[[i]] <- matrix(Q[[i]][, 1:d[i]], p, d[i])
  
  
  #Calculation of the remaining parameters of the selected model
  
  # PARAMETER a
  ai <- matrix(NA, K, max(d))
  if ( model%in%c('AKJBKQKDK', 'AKJBQKDK', 'AKJBKQKD', 'AKJBQKD') ){
    for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
  } else if ( model%in%c('AKBKQKDK', 'AKBQKDK' , 'AKBKQKD', 'AKBQKD') ){
    for (i in 1:K) ai[i, ] <- rep(sum(ev[i, 1:d[i]])/d[i], length=max(d))
  } else if(model=="AJBQD"){
    for (i in 1:K) ai[i, ] <- ev[1:d[1]]
  } else if(model=="ABQD") {
    ai[] <- sum(ev[1:d[1]])/d[1]
  } else {
    a <- 0
    eps <- sum(prop*d)
    for (i in 1:K) a <- a + sum(ev[i, 1:d[i]])*prop[i]
    ai <- matrix(a/eps, K, max(d))
  }
  
  # PARAMETER b
  bi <- c()
  denom = min(N,p)
  if ( model%in%c('AKJBKQKDK', 'AKBKQKDK', 'ABKQKDK', 'AKJBKQKD', 'AKBKQKD', 'ABKQKD') ){
    for(i in 1:K){
      remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
      # bi[i] <- sum(ev[i, (d[i]+1):min(N, p)])/(p-d[i])
      bi[i] <- remainEV/(p-d[i]) #pour moi c'est p au lieu de denom
    }
  } else if ( model%in%c("ABQD", "AJBQD") ){
    remainEV = traceVect - sum(ev[1:d[1]])
    # bi[1:K] <- sum(ev[(d[1]+1):min(N, p)])/(min(N, p)-d[1])
    bi[1:K] <- remainEV/(denom-d[1])
  } else {
    b <- 0
    eps <- sum(prop*d)
    for(i in 1:K){
      remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
      # b <- b + sum(ev[i, (d[i]+1):min(N, p)])*prop[i]
      b <- b + remainEV*prop[i]
    }
    bi[1:K] <- b/(min(N, p)-eps)
  }
  
  # We adjust the values of b if they are too low
  #bi[bi<noise.ctrl] = noise.ctrl
  #---------------------------------------
  
  #A@5 new add formula 17 from the paper in here: These are the alpha's if the numbers are not larger than 0.5 must use optimize
  
  
  etamax1 <- function(eta, zk, vk, x, muk, wk, Qki, aki, bki){
    p=ncol(x)
    N=nrow(x)
    vki = matrix(vk, N, 1, byrow=TRUE)
    zki = matrix(zk, N, 1, byrow=TRUE) 
    muki = matrix(muk, 1, p, byrow=TRUE)
    cone=.imahalanobis(x, muki, wk, Qki, aki, bki)
    return(-p/2*sum(zki*(1-vki)*log(eta)) - 1/(2*eta)*sum(zki*(1-vki)*cone))
  }# A@5 what value should tol have, we are using from the original paper
  for(i in 1:K){
    exv = matrix(vx[,i], N, 1, byrow=TRUE)
    exz = matrix(t[,i], N, 1, byrow=TRUE) 
    exQki = Q[[i]]
    exak = sqrt(diag(1/ai[i,1:d[i]],d[i]))
    exbk = bi[i]
    exmuk = matrix(mu[i,], 1, p, byrow=TRUE)
    exwk = fpcaobj[[i]]$W
    excone=.imahalanobis(x, exmuk, exwk, exQki, exak, exbk)
    #print(paste0("imahSum ",sum(exz*(1.0-exv)*excone) ))
    #print(paste0("No imahSum ",sum(exz*(1.0-exv)) )
    ala<-sum(exz*(1.0-exv))
    blb<-sum(exz*(1.0-exv)*excone)
    etry<-blb/(p*ala)
    etamax=1000
    #if ((etry>1)&&(etry<1000)) etax[i]=etry
    #else {if(etry<1) etax[i]<-1.01 else etax[i]<-100}
    #    print(paste0("caletax=",etax[i] ))
    ##############################
    # if(etry < 1) {
    ## etax[i] <- optimize(etamax1, c(1,100), maximum = TRUE,zk = exz,vk = exv,x = x, muk = exmuk,wk = exwk,Qki =exQki , aki =exak , bki =exbk)$maximum
    #  etax[i]=1.01
    #  } else if(etry <etamax)
    #   etax[i]=etry
    ################################
    if(etry < 1)
      etax[i]=1.01
    else
    { if((etry > etamax)||(is.nan(etry))) {
      etax[i] <- optimize(etamax1, c(1,etamax), maximum = TRUE,zk = exz,vk = exv,x = x, muk = exmuk,wk = exwk,Qki =exQki , aki =exak , bki =exbk)$maximum
      
    } else 
      etax[i]=etry
    }}
  #A@5 add the CM2 step
  
  # A@5 new mahalanobis rowSums(A^2)+1/b[i]*rowSums(B^2) from e step
  ###########new mu in terms of new eta ca in CNMixt
  # correctionX <- array(1,c(N,K),dimnames=list(1:N,paste("comp.",1:K,sep=""))) 
  
  #correctionX <- vx+(1.0-vx)*matrix(1.0/etax,N,K,byrow=TRUE)#like fact in CNMixt  
  #mu <- matrix(NA, K, p)
  #corX=matrix(0,N,K)
  #corX=t*correctionX
  #for (i in 1:K) {
  #mu[i, ] <- colSums(x*t[, i])/n[i] #A@5 add the big parenthesis from formula 18 and becomes something like colSums(x*t[,i]* a slive of v)
  #A@5 and n[i] must be multiplied by some form of colsum of the product between t and the same thing as above
  # a@5 new
  
  # mu[i,]     <- colSums(x*corX[,i])/sum(corX[,i])
  #}
  ###################################
  #A@5 check the other library on how to use optimize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  list(model=model, K=K, d=d, a=ai, b=bi, mu=mu, prop=prop, alphax=alphax,etax=etax,ev=ev, Q=Q,fpcaobj = fpcaobj) #A@5 add the extra params(alpha and eta)
  
  
}

.mypca.fd1 <- function(fdobj, Ti,CorI){
  if (class(fdobj)=="list"){ #A@5 NOT OUR CASE 
    #sauvegarde de la moyenne avant centrage
    mean_fd<-list()
    for (i in 1:length(fdobj)){
      mean_fd[[i]]<-fdobj[[i]]
    }
    
    #centrage des objets fonctionnels
    for (i in 1:length(fdobj)){
      coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj[[i]]$coefs))) * fdobj[[i]]$coefs, 1, sum) / sum(Ti) 
      fdobj[[i]]$coefs <- sweep(fdobj[[i]]$coefs, 1, coefmean)  
      mean_fd[[i]]$coefs = as.matrix(data.frame(mean=coefmean)) 
    }
    
    #Cr?ation des matrices des produits des bases de fonction 
    
    for (i in 1:length(fdobj)){
      name<-paste('W_var',i,sep='')
      W_fdobj<-inprod(fdobj[[i]]$basis,fdobj[[i]]$basis)
      assign(name,W_fdobj)
    }
    
    #Ajout des 0 ? gauche et ? droite des matrices W avant leur fusion en matrice phi
    prow<-dim(W_fdobj)[[1]]
    pcol<-length(fdobj)*prow
    W1<-cbind(W_fdobj,matrix(0,nrow=prow,ncol=(pcol-ncol(W_fdobj))))
    W_list<-list()
    for (i in 2:(length(fdobj))){
      W2<-cbind(matrix(0,nrow=prow,ncol=(i-1)*ncol(W_fdobj)),get(paste('W_var',i,sep='')),matrix(0,nrow=prow,ncol=(pcol-i*ncol(W_fdobj))))
      W_list[[i-1]]<-W2
    }
    
    #Cr?ation de la matrice phi
    W_tot<-rbind(W1,W_list[[1]])
    if (length(fdobj)>2){
      for(i in 2:(length(fdobj)-1)){
        W_tot<-rbind(W_tot,W_list[[i]])
      }
    }
    W_tot[W_tot<1e-15]=0
    
    
    #Cr?ation de la matrice de coef
    coef<-t(fdobj[[1]]$coefs)
    for (i in 2:length(fdobj)){
      coef<-cbind(coef,t(fdobj[[i]]$coefs))
    }
    
    #mat_interm<-1/sqrt((sum(Ti)-1))*coef%*%(W_tot)^(1/2)
    
    #D?but Partie mise en commentaire
    #mat_interm<-1/sqrt((sum(Ti)))*coef%*%chol(W_tot,pivot=TRUE)
    #cov<-(.repmat(Ti,n=pcol,p=1)*t(mat_interm))%*%mat_interm
    #Fin de partie mise en commentaire
    
    #D?but correction
    #Construction matrice triangulaire de Choleski
    W_m <-  chol(W_tot) #A@5 this is sqrt(W)
    #Matrice de covariance 
    mat_cov <- crossprod(t(.repmat(sqrt(Ti),n=dim(t(coef))[[1]],p=1)*t(coef)))/sum(Ti) 
    cov = W_m %*% mat_cov %*% t(W_m) 
    #Fin correction
    
    
    valeurs<-Eigen(cov)
    valeurs_propres<-valeurs$values
    vecteurs_propres<-valeurs$vectors
    #bj<-solve((W_tot)^(1/2))%*%vecteurs_propres
    bj<-solve(chol(W_tot))%*%vecteurs_propres
    fonctionspropres<-fdobj[[1]]
    fonctionspropres$coefs<-bj
    scores<-coef%*%W_tot%*%bj
    
    varprop<-valeurs_propres/sum(valeurs_propres)
    
    pcafd<-list(valeurs_propres=valeurs_propres,harmonic=fonctionspropres,scores=scores,covariance=cov,U=bj,varprop=varprop,meanfd=mean_fd,W=W_tot)
    
    
  }else if (class(fdobj)!="list") { #A@5 THSI IS OUR CASE 
    #Calcul de la moyenne par groupe
    mean_fd<-fdobj
    #Centrer les objets fonctionnels par groupe
    
    coefmean <- apply(t(as.matrix(CorI) %*% matrix(1,1,nrow(fdobj$coefs))) * fdobj$coefs, 1, sum) / sum(CorI) #A@5 here they compute miu(i) - second line of fomula 19
    # coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj$coefs))) * fdobj$coefs, 1, sum) / sum(Ti) #A@5 here they compute miu(i) - second line of fomula 19
    
    fdobj$coefs <- sweep(fdobj$coefs, 1, coefmean) #A@5 what does sweep do?substracts coefmean from gamma, i.e gammai-mui
    mean_fd$coefs = as.matrix(data.frame(mean=coefmean))
    
    #Calcul de la matrice des produits scalaires #A@5 this stays the SAME CAN"T TOUCH THIS!
    #----------------------------------------------------------
    W<-inprod(fdobj$basis,fdobj$basis)
    #pour ?viter les soucis num?riques on arrondit ? 0 les produits scalaires tr?s petits
    W[W<1e-15]=0
    #on centre les coefficients
    coef<-t(fdobj$coefs)
    #mat_interm<-1/sqrt((sum(Ti)-1))*coef%*%(W)^(1/2)
    
    #D?but partie mise en commentaire
    #mat_interm<-1/sqrt((sum(Ti)))*coef%*%chol(W)
    #cov<-(.repmat(Ti,n=dim(W)[[1]],p=1)*t(mat_interm))%*%mat_interm
    #Fin de la partie mise en commentaire
    
    #D?but correction
    #Construction matrice triangulaire de Choleski
    W_m <-  chol(W)
    #-----------------------------------------------------------   
    #A@5    
    
    #Matrice de covariance #A@5 mat_cov is formula 19 without the big bracket which we have to add
    mat_cov <- crossprod(t(.repmat(sqrt(CorI),n=dim(t(coef))[[1]],p=1)*t(coef)))/sum(Ti) 
    cov = W_m %*% mat_cov %*% t(W_m) #A@5 this is the last formula on page 7
    #Fin correction
    #A@5 CAN'T TOUCH THIS
    #--------------------------------------------  
    valeurs<-Eigen(cov)
    valeurs_propres<-valeurs$values
    vecteurs_propres<-valeurs$vectors
    #Calcul de U
    #U<-chol(W)%*%vecteurs_propres
    #Calcul des coefficients des fonctions propres
    fonctionspropres<-fdobj
    bj<-solve(chol(W))%*%vecteurs_propres
    fonctionspropres$coefs<-bj
    #calcul des scores selon la formule de pca.fd
    scores<-inprod(fdobj,fonctionspropres)
    
    varprop<-valeurs_propres/sum(valeurs_propres)
    
    pcafd <-list(valeurs_propres=valeurs_propres,harmonic=fonctionspropres,scores=scores,covariance=cov,U=bj,meanfd=mean_fd,W=W)
    #-------------------------------------------    
  }
  class(pcafd) <- "pca.fd"
  return(pcafd)
}


.hddc_ari <- function(x,y){
  #This function is drawn from the mclust package
  x <- as.vector(x)
  y <- as.vector(y)
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}



####
#### CONTROLS ####
####

.hddc_control = function(call){
  
  prefix = "HDDC: "
  .myCallAlerts(call, "data", "list,fd", 3, TRUE, prefix)
  .myCallAlerts(call, "K", "integerVector", 3, FALSE, prefix)
  .myCallAlerts(call, "model", "vector", 3, FALSE, prefix)
  .myCallAlerts(call, "threshold", "numericVectorGE0LE1", 3, FALSE, prefix)
  .myCallAlerts(call, "criterion", "character", 3, FALSE, prefix)
  .myCallAlerts(call, "com_dim", "singleIntegerGE1", 3, FALSE, prefix)
  .myCallAlerts(call, "itermax", "singleIntegerGE0", 3, FALSE, prefix)
  .myCallAlerts(call, "eps", "singleNumericGE0", 3, FALSE, prefix)
  .myCallAlerts(call, "graph", "singleLogical", 3, FALSE, prefix)
  .myCallAlerts(call, "algo", "singleCharacter", 3, FALSE, prefix)
  .myCallAlerts(call, "d_select", "singleCharacter", 3, FALSE, prefix)
  .myCallAlerts(call, "init", "singleCharacter", 3, FALSE, prefix)
  .myCallAlerts(call, "show", "singleLogical", 3, FALSE, prefix)
  .myCallAlerts(call, "mini.nb", "integerVectorGE1", 3, FALSE, prefix)
  .myCallAlerts(call, "scaling", "singleLogical", 3, FALSE, prefix)
  .myCallAlerts(call, "min.individuals", "singleIntegerGE2", 3, FALSE, prefix)
  .myCallAlerts(call, "noise.ctrl", "singleNumericGE0", 3, FALSE, prefix)
  .myCallAlerts(call, "mc.cores", "singleIntegerGE1", 3, FALSE, prefix)
  .myCallAlerts(call, "nb.rep", "singleIntegerGE1", 3, FALSE, prefix)
  .myCallAlerts(call, "keepAllRes", "singleLogical", 3, FALSE, prefix)
  .myCallAlerts(call, "d_max", "singleIntegerGE1", 3, FALSE, prefix)
  
  
  ####
  #### SPECIFIC controls
  ####
  
  # Getting some elements
  data = eval.parent(call[["data"]], 2)
  K = eval.parent(call[["K"]], 2)
  init = eval.parent(call[["init"]], 2)
  criterion = eval.parent(call[["criterion"]], 2)
  
  if (is.fd(data)){
    # No NA in the data:
    if (any(is.na(data$coefs))) stop("NA values in the data are not supported. Please remove them beforehand.")
    # Size of the data
    if(any(K>2*NROW(data$coefs))) stop("The number of observations must be at least twice the number of clusters ")
    
  }else{
    # No NA in the data:
    for (i in 1:length(data)){
      if (any(is.na(data[[i]]$coefs))) stop("NA values in the data are not supported. Please remove them beforehand.")
    }
    # Size of the data
    if(any(K>2*NROW(data[[1]]$coefs))) stop("The number of observations must be at least twice the number of clusters ")
  }
  
  # Initialization Controls
  if(!is.null(init)){
    
    # we get the value of the initialization
    init = .myAlerts(init, "init", "singleCharacterMatch.arg", "HDDC: ", c('random', 'kmeans', 'mini-em', 'param', "vector"))
    
    # Custom initialization => controls and setup
    if(init == "vector"){
      fdobj = data
      if (class(fdobj)!='list') {x = t(fdobj$coefs)}
      else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
      .myCallAlerts(call, "init.vector", "(integer,factor)Vector", 3, FALSE, prefix)
      
      init.vector = eval.parent(call[["init.vector"]], 2)
      
      if(is.null(init.vector)) stop("HDDC: When init='vector', the argument 'init.vector' should be provided.")
      
      if(length(unique(K))>1) stop("HDDC: Several number of classes K cannot be estimated when init='vector'.")
      
      init.vector <- unclass(init.vector)
      if(K!=max(init.vector)) stop("The number of class K, and the number of classes in the initialization vector are different")
      
      if( length(init.vector)!=nrow(x) ) stop("The size of the initialization vector is different of the size of the data")
    }
    
    # The param init
    if (init=='param' && nrow(data)<ncol(data)){
      stop("The 'param' initialization can't be done when N<p")
    }
    
    # The mini.em init
    if (init=='mini-em'){
      
      mini.nb = eval.parent(call[["mini.nb"]], 2)
      
      if(!is.null(mini.nb) && length(mini.nb)!=2){
        stop("The parameter mini.nb must be a vector of length 2 with integers\n")
      }
      
    }
  }
  
}

.default_kmeans_control = function(control){
  
  .myAlerts(control,"kmeans.control","list","kmeans controls: ")
  
  #
  # Default values of the control parameters
  #
  
  myDefault = list()
  myDefault$iter.max = 10
  myDefault$nstart = 1
  myDefault$algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen")
  myDefault$trace = FALSE
  
  #
  # Types of each arg
  #
  
  myTypes = c("singleIntegerGE1", "singleIntegerGE1", "match.arg", "singleLogical")
  
  #
  # Recreation of the kmeans controls + Alerts
  #
  
  control = .matchTypeAndSetDefault(control, myDefault, myTypes, "kmeans list of controls: ")
  
  return(control)
}

#=================================#
# This file contains all the
# "control" functions
#=================================#

# Possible elements of myAlerts:
#

.myCallAlerts = function(call, name, myType, nParents=1, mustBeThere=FALSE, prefix=""){
  # This function basically calls the function myAlerts, but the arguments are different
  
  if( name %in% names(call) ){
    # we check the element exists => to provide a fine error
    what = call[[name]]
    val = try(eval.parent(what, nParents), silent = TRUE)
    # browser()
    if( "try-error" %in% class(val) ){
      if( class(what)=="name" ){
        # it means the variable was not found
        stop(prefix,"For argument '",name,"': object '",what,"' not found.", call. = FALSE)
      } else {
        stop(prefix,"For argument '",name,"': expression ",as.character(as.expression(what))," could not be evaluated.", call. = FALSE)
      }
      
    } else {
      a = .myAlerts(val, name, myType, prefix)
      return(a)
    }
  } else if(mustBeThere) {
    stop(prefix, "The argument '", name, "' must be provided.", call. = FALSE)
  }
}

.myAlerts = function(x, name, myType, prefix="", charVec){
  # Format of my types:
  #   - single => must be of lenght one
  #   - Vector => must be a vector
  #   - Matrix => must be a matrix
  #   - GE/GT/LE/LT: greater/lower than a given value
  #   - predefinedType => eg: numeric, integer, etc
  #   - match.arg => very specific => should match the charVec
  # If there is a parenthesis => the class must be of specified types:
  # ex: "(list, data.frame)" must be a list of a data.frame
  
  ignore.case = TRUE
  
  firstMsg = paste0(prefix,"The argument '",name,"' ")
  
  # simple function to extract a pattern
  # ex: if my type is VectorIntegerGE1 => myExtract("GE[[:digit:]]+","VectorIntegerGE1") => 1
  myExtract = function(expr, text, trim=2){
    start = gregexpr(expr,text)[[1]] + trim
    length = attr(start,"match.length") - trim
    res = substr(text,start,start+length-1)
    as.numeric(res)
  }
  
  #
  # General types handling
  #
  
  loType = tolower(myType)
  
  if(grepl("single",loType)){
    if(length(x)!=1) stop(firstMsg,"must be of length one.", call. = FALSE)
  }
  
  if(grepl("vector",loType) && !grepl("factor",loType)){
    if(!is.vector(x)) stop(firstMsg,"must be a vector.", call. = FALSE)
    if(is.list(x)) stop(firstMsg,"must be a vector (and not a list).", call. = FALSE)
  }
  
  res = .checkTheTypes(loType, x)
  if(!res$OK) stop(firstMsg,res$message, call. = FALSE)
  
  # # INTEGER is a restrictive type that deserves some explanations (not included in getTheTypes)
  # if(grepl("integer",loType)){
  #     if(grepl("single",loType)){
  #         if(!is.numeric(x)) stop(firstMsg,"must be an integer (right now it is not even numeric).", call. = FALSE)
  #         if(!(is.integer(x) || x%%1==0)) stop(firstMsg,"must be an integer.", call. = FALSE)
  #     } else {
  #         if(!is.numeric(x)) stop(firstMsg,"must be composed of integers (right now it is not even numeric).", call. = FALSE)
  #         if(!(is.integer(x) || all(x%%1==0))) stop(firstMsg,"must be composed of integers.", call. = FALSE)
  #     }
  # }
  
  # GE: greater or equal // GT: greater than // LE: lower or equal // LT: lower than
  if(grepl("ge[[:digit:]]+",loType)){
    n = myExtract("ge[[:digit:]]+", loType)
    if( !all(x>=n) ) stop(firstMsg,"must be greater than, or equal to, ", n,".", call. = FALSE)
  }
  if(grepl("gt[[:digit:]]+",loType)){
    n = myExtract("gt[[:digit:]]+", loType)
    if( !all(x>n) ) stop(firstMsg,"must be strictly greater than ", n,".", call. = FALSE)
  }
  if(grepl("le[[:digit:]]+",loType)){
    n = myExtract("le[[:digit:]]+", loType)
    if( !all(x<=n) ) stop(firstMsg,"must be lower than, or equal to, ",n,".", call. = FALSE)
  }
  if(grepl("lt[[:digit:]]+",loType)){
    n = myExtract("lt[[:digit:]]+", loType)
    if( !all(x<n) ) stop(firstMsg,"must be strictly lower than ", n,".", call. = FALSE)
  }
  
  #
  # Specific Types Handling
  #
  
  if(grepl("match.arg",loType)){
    if(ignore.case){
      x = toupper(x)
      newCharVec = toupper(charVec)
    } else {
      newCharVec = charVec
    }
    
    if( is.na(pmatch(x, newCharVec)) ){
      n = length(charVec)
      if(n == 1){
        msg = paste0("'",charVec,"'")
      } else {
        msg = paste0("'", paste0(charVec[1:(n-1)], collapse="', '"), "' or '",charVec[n],"'")
      }
      stop(firstMsg, "must be one of:\n", msg, ".", call. = FALSE)
    } else {
      qui = pmatch(x, newCharVec)
      return(charVec[qui])
    }
  }
}

.matchTypeAndSetDefault = function(myList, myDefault, myTypes, prefix){
  # Cette fonction:
  #   i) verifie que tous les elements de la liste sont valides
  #   ii) mes les valeurs par defauts si elles certaines valeurs sont manquantes
  #   iii) Envoie des messages d'erreur si les typages ne sont pas bons
  # En fait cette fonction "coerce" myList en ce qu'il faudrait etre (donne par myDefault)
  
  # 1) check that the names of the list are valid
  if(is.null(myList)) myList = list()
  list_names = names(myList)
  
  if(length(list_names)!=length(myList) || any(list_names=="")){
    stop(prefix,"The elements of the list should be named.", call. = FALSE)
  }
  
  obj_names = names(myDefault)
  
  isHere = pmatch(list_names,obj_names)
  
  if(anyNA(isHere)){
    if(sum(is.na(isHere))==1) stop(prefix, "The following argument is not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
    else stop(prefix, "The following arguments are not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
  }
  
  # 2) We set the default values and run Warnings
  res = list()
  for(i in 1:length(obj_names)){
    obj = obj_names[i]
    qui = which(isHere==i) # qui vaut le numero de l'objet dans myList
    type = myTypes[i] # we extract the type => to control for "match.arg" type
    if(length(qui)==0){
      # we set to the default if it's missing
      if(type == "match.arg") {
        res[[obj]] = myDefault[[i]][1]
      } else {
        res[[obj]] = myDefault[[i]]
      }
    } else {
      # we check the integrity of the value
      val = myList[[qui]]
      if(type == "match.arg"){
        # If the value is to be a match.arg => we use our controls and not
        # directly the one of the function match.arg()
        charVec = myDefault[[i]]
        .myAlerts(val, obj, "singleCharacterMatch.arg", prefix, charVec)
        val = match.arg(val, charVec)
      } else {
        .myAlerts(val, obj, type, prefix)
      }
      
      res[[obj]] = val
    }
  }
  
  return(res)
}



.checkTheTypes = function(str, x){
  # This function takes in a character string describing the types of the
  # element x => it can be of several types
  
  # types that are controlled for:
  allTypes = c("numeric", "integer", "character", "logical", "list", "data.frame", "matrix", "factor")
  
  OK = FALSE
  message = c()
  
  for(type in allTypes){
    
    if(grepl(type, str)){
      
      # we add the type of the control
      message = c(message, type)
      
      if(type == "numeric"){
        if(!OK & is.numeric(x)){
          OK = TRUE
        }
      } else if(type == "integer"){
        if(is.numeric(x) && (is.integer(x) || all(x%%1==0))){
          OK = TRUE
        }
      } else if(type == "character"){
        if(is.character(x)){
          OK = TRUE
        }
      } else if(type == "logical"){
        if(is.logical(x)){
          OK = TRUE
        }
      } else if(type == "list"){
        if(is.list(x)){
          OK = TRUE
        }
      } else if(type == "data.frame"){
        if(is.data.frame(x)){
          OK=TRUE
        }
      } else if(type == "matrix"){
        if(is.matrix(x)){
          OK = TRUE
        }
      } else if(type == "factor"){
        if(is.factor(x)){
          OK = TRUE
        }
      }
    }
    
    if(OK) break
  }
  
  if(length(message) == 0) OK = TRUE #ie there is no type to be searched
  else if(length(message) >= 3){
    n = length(message)
    message = paste0("must be of type: ",  paste0(message[1:(n-1)], collapse = ", "), " or ", message[n], ".")
  } else {
    message = paste0("must be of type: ",  paste0(message, collapse = " or "), ".")
  }
  
  
  return(list(OK=OK, message=message))
}




.hdclassif_dim_choice <- function(ev, n, method, threshold, graph, noise.ctrl){
  # Selection of the intrinsic dimension
  # browser()
  N <- sum(n)
  prop <- n/N
  K = ifelse(is.matrix(ev), nrow(ev), 1)
  
  # browser()
  
  if(is.matrix(ev) && K>1){
    p <- ncol(ev)
    if(method=="cattell"){
      dev <- abs(apply(ev, 1, diff))
      max_dev <- apply(dev, 2, max, na.rm=TRUE)
      dev <- dev/rep(max_dev, each=p-1)
      d <- apply((dev>threshold)*(1:(p-1))*t(ev[, -1]>noise.ctrl), 2, which.max)
      
      if(graph){
        op = par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9), 1+floor(K/4)-1*(K==12)+1*(K==7)))
        for(i in 1:K){
          sub1 <- paste("Class #", i, ",  d", i, "=", d[i], sep="")
          Nmax <- max(which(ev[i, ]>noise.ctrl))-1
          plot(dev[1:(min(d[i]+5, Nmax)), i], type="l", col="blue", main=paste("Cattell's Scree-Test\n", sub1, sep=""), ylab=paste("threshold =", threshold), xlab="Dimension", ylim=c(0, 1.05))
          abline(h=threshold, lty=3)
          points(d[i], dev[d[i], i], col='red')
        }
        par(op)
      }
    } else if(method=="bic"){
      
      d <- rep(0, K)
      if(graph) op = par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9), 1*(1+floor(K/4)-1*(K==12)+1*(K==7))))
      
      for (i in 1:K) {
        B <- c()
        Nmax <- max(which(ev[i, ]>noise.ctrl))-1
        p2 <- sum(!is.na(ev[i, ]))
        Bmax <- -Inf
        for (kdim in 1:Nmax){
          if ((d[i]!=0 & kdim>d[i]+10)) break
          a <- sum(ev[i, 1:kdim])/kdim
          b <- sum(ev[i, (kdim+1):p2])/(p2-kdim)
          if (b<0 | a<0){
            B[kdim] <- -Inf
          } else {
            L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
            B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
          }
          
          if ( B[kdim]>Bmax ){
            Bmax <- B[kdim]
            d[i] <- kdim
          }
        }
        
        if(graph){
          plot(B, type='l', col=4, main=paste("class #", i, ",  d=", d[i], sep=''), ylab='BIC', xlab="Dimension")
          points(d[i], B[d[i]], col=2)
        }
      }
      if(graph) par(op)
    }
  } else{
    ev <- as.vector(ev)
    p <- length(ev)
    
    if(method=="cattell"){
      dvp <- abs(diff(ev))
      Nmax <- max(which(ev>noise.ctrl))-1
      if (p==2) d <- 1
      else d <- max(which(dvp[1:Nmax]>=threshold*max(dvp[1:Nmax])))
      diff_max <- max(dvp[1:Nmax])
      
      if(graph){
        plot(dvp[1:(min(d+5, p-1))]/diff_max, type="l", col="blue", main=paste("Cattell's Scree-Test\nd=", d, sep=''), ylab=paste("threshold =", threshold, sep=' '), xlab='Dimension', ylim=c(0, 1.05))
        abline(h=threshold, lty=3)
        points(d, dvp[d]/diff_max, col='red')
      }
    } else if(method=="bic"){
      d <- 0
      Nmax <- max(which(ev>noise.ctrl))-1
      B <- c()
      Bmax <- -Inf
      for (kdim in 1:Nmax){
        if (d!=0 && kdim>d+10) break
        a <- sum(ev[1:kdim])/kdim
        b <- sum(ev[(kdim+1):p])/(p-kdim)
        if (b<=0 | a<=0) B[kdim] <- -Inf
        else{
          L2 <- -1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))*N
          B[kdim] <- 2*L2 - (p+kdim*(p-(kdim+1)/2)+1)*log(N)
        }
        if ( B[kdim]>Bmax ){
          Bmax <- B[kdim]
          d <- kdim
        }
      }
      
      if(graph){
        plot(B, type='l', col=4, main=paste("BIC criterion\nd=", d, sep=''), ylab='BIC', xlab="Dimension")
        points(d, B[d], col=2)
      }
    }
  }
  return(d)
}

.hdclassif_bic <- function(par, p, data=NULL){
  model <- par$model
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  N <- par$N
  prop <- par$prop
  
  if(length(b)==1){
    #update of b to set it as variable dimension models
    eps <- sum(prop*d)
    n_max <- if(model%in%c("ABQD", "AJBQD")) length(par$ev) else ncol(par$ev)
    b <- b*(n_max-eps)/(p-eps)
    b <- rep(b, length=K)
  }
  if (length(a)==1) a <- matrix(a, K, max(d))
  else if (length(a)==K) a <- matrix(a, K, max(d))
  else if (model=='AJBQD') a <- matrix(a, K, d[1], byrow=TRUE)
  
  if(min(a, na.rm=TRUE)<=0 | any(b<0)) return(-Inf)
  
  if (is.null(par$loglik)){
    som_a <- c()
    #probabil likelihoodul initial
    for (i in 1:K) som_a[i] <- sum(log(a[i, 1:d[i]]))
    L <- -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop)+ p*(1+log(2*pi))))*N
  }
  #@ we don't care about these models
  else if (model%in%c("ABQD", "AJBQD")){
    Q <- rep(list(par$Q), K)
    K_pen <- matrix(0, K, N)
    for (i in 1:K) {
      s <- sum(log(a[i, 1:d[i]]))
      X <- data-matrix(mu[i, ], N, p, byrow=TRUE)
      proj <- (X%*%Q[[i]])%*%t(Q[[i]])
      A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i, 1:d[i]], d[i]))
      B <- X-proj
      K_pen[i, ] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
    }
    A <- -1/2*t(K_pen)
    L <- sum(log(rowSums(exp(A-apply(A, 1, max))))+apply(A, 1, max))
  } else L <- par$loglik[length(par$loglik)]
  
  #@ add 2K parameters for alpha and eta in ro
  ro <- K*(p+2)+K-1
  tot <- sum(d*(p-(d+1)/2))
  D <- sum(d)
  d <- d[1]
  to <- d*(p-(d+1)/2)
  if (model=='AKJBKQKDK') m <- ro+tot+D+K
  else if (model=='AKBKQKDK') m <- ro+tot+2*K
  else if (model=='ABKQKDK') m <- ro+tot+K+1
  else if (model=='AKJBQKDK') m <- ro+tot+D+1
  else if (model=='AKBQKDK') m <- ro+tot+K+1
  else if (model=='ABQKDK') m <- ro+tot+2
  bic <- -(-2*L+m*log(N))
  
  #calcul ICL
  t = par$posterior
  #if(!is.null(t)){
  # means we are in HDDC
  Z = ((t - apply(t, 1, max))==0) + 0
  icl = bic - 2*sum(Z*log(t+1e-15))
  # } else {
  #   # Si HDDA, entropie est nulle => car classes pures
  #   icl = bic
  #}
  
  return(list(bic = bic, icl = icl))
}

.hdc_getComplexity = function(par, p){
  model <- par$model
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  N <- par$N
  prop <- par$prop
  #@ add in ro 2K parameters for alpha and eta
  ro <- K*(p+2)+K-1
  tot <- sum(d*(p-(d+1)/2))
  D <- sum(d)
  d <- d[1]
  to <- d*(p-(d+1)/2)
  if (model=='AKJBKQKDK') m <- ro+tot+D+K
  else if (model=='AKBKQKDK') m <- ro+tot+2*K
  else if (model=='ABKQKDK') m <- ro+tot+K+1
  else if (model=='AKJBQKDK') m <- ro+tot+D+1
  else if (model=='AKBQKDK') m <- ro+tot+K+1
  else if (model=='ABQKDK') m <- ro+tot+2
  
  return(m)
}

.hdc_getTheModel = function(model, all2models = FALSE){
  # Function used to get the models from number or names
  
  model_in = model
  
  if(!is.vector(model)) stop("The argument 'model' must be a vector.")
  
  if(anyNA(model)) stop("The argument 'model' must not contain any NA.")
  
  ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD")
  
  model = toupper(model)
  
  if(length(model)==1 && model=="ALL"){
    if(all2models) model <- 1:14
    else return("ALL")
  }
  
  qui = which(model %in% 1:14)
  model[qui] = ModelNames[as.numeric(model[qui])]
  
  # We order the models properly
  qui = which(!model%in%ModelNames)
  if (length(qui)>0){
    if(length(qui)==1){
      msg = paste0("(e.g. ", model_in[qui], " is incorrect.)")
    } else {
      msg = paste0("(e.g. ", paste0(model_in[qui[1:2]], collapse=", or "), " are incorrect.)")
    }
    stop("Invalid model name ", msg)
  }
  
  # warning:
  if(max(table(model))>1) warning("The model vector, argument 'model', is made unique (repeated values are not tolerated).")
  
  mod_num <- c()
  for(i in 1:length(model)) mod_num[i] <- which(model[i]==ModelNames)
  mod_num <- sort(unique(mod_num))
  model <- ModelNames[mod_num]
  
  return(model)
}



####
#### Utilities ####
####


.addCommas = function(x) sapply(x, .addCommas_single )

.addCommas_single = function(x){
  # Cette fonction ajoute des virgules pour plus de
  # visibilite pour les (tres longues) valeurs de vraisemblance
  
  if(!is.finite(x)) return(as.character(x))
  
  s = sign(x)
  x = abs(x)
  
  decimal = x - floor(x)
  if(decimal>0) dec_string = substr(decimal, 2, 4)
  else dec_string = ""
  
  entier = as.character(floor(x))
  
  quoi = rev(strsplit(entier, "")[[1]])
  n = length(quoi)
  sol = c()
  for(i in 1:n){
    sol = c(sol, quoi[i])
    if(i%%3 == 0 && i!=n) sol = c(sol, ",")
  }
  
  res = paste0(ifelse(s==-1, "-", ""), paste0(rev(sol), collapse=""), dec_string)
  res
}


.repmat <- function(v,n,p){ #A@5 WHAT IS THIS???????????????????????
  if (p==1){M = cbind(rep(1,n)) %*% v} #A@5 a matrix of column of v
  else { cat('!'); M = matrix(rep(v,n),n,(length(v)*p),byrow=T)}
  M
}

.diago <- function(v){
  if (length(v)==1){ res = v }
  else { res = diag(v)}
  res
}
#imahalanobis(x, exmuk, exwk, exQki, exak, exbk)
.imahalanobis <- function(x, muk, wk, Qk, aki,bki) {
  # w = fpcaobj[[i]]$W
  # *************************************************************************** #
  # should be called as imahalanobis(x, mu[i,], w[i,], Q[[i]], a[i,], d[i], b[i])
  # return formula on top of page 9
  # *************************************************************************** #
  
  # X <- x - matrix(mu[i,], N, p, byrow=TRUE)
  #  Qi = par$fpcaobj[[i]]$W %*% Q[[i]]
  # proj <- (X%*%Qi)%*%t(Qi)
  #  A <- (-proj)%*%Qi%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
  # B <- X-proj
  #K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
  
  # Code modified to take vectors instead of matrices
  p <- ncol(x)
  N <- nrow(x)
  res<-rep(0,N)
  #X <- x - matrix(muk[i,], N, p, byrow=TRUE)
  X <- x - matrix(muk, N, p, byrow=TRUE)
  #Qi <- wk %*% Qk[[i]]
  Qi <- wk %*% Qk
  #  Qi <- Qk
  proj <- (X%*%Qi)%*%t(Qi)
  # A <- (-proj)%*%Qi%*%sqrt(diag(1/ak[i,1:dk[i]],dk[i])) #A@5 NORM of this is the first term in formula 15
  A <- (-proj)%*%Qi%*%aki
  B <- X-proj #A@5 NOMR of this is the second term in formula 15
  # sum of the first two terms is new mahalanobis
  #return(rowSums(A^2)+1/bk[i]*rowSums(B^2))
  res<-rowSums(A^2)+rowSums(B^2)/bki
  return(res)
  
}
#ianos<-cluster_comparison(1,K=3, fun_names=c("cfunHDDC"), nrep=1, eta=c(5, 50, 15), nsplines=35, ncurves=999, alpha=c(0.9, 0.9, 0.9),threshold=0.2, alphamin=c(0.75))

#plot(data$fd[1:333], col=CL[data$groupd[1:333]])
#plot(data$fd[667:999], col=CL[data$groupd[667:999]])
#plot(data$fd[334:666], col=CL[data$groupd[334:666]])
#ianos1=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))
#83/115
#ianos2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.75))
#86/115
#ianos3=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=100,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.5))#78/115
#ianos=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=100,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#80/115
#ianos2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("funHDDC"), nrep=100,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#74/115
#ianos4=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#76/115
#ianos5=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.5))#77/115
#ianos6=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.75))#76/115
#ianos7=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.65))#77/115
#ianos8=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.9))#76/115
#ianos9=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.95))#77/115
#in1=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#80/115 16 outliers
#in2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))#80/115 13 outliers
#in3=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.01,alphamin=c(0.85))#79/115 4 outliers
#in4=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.85))#69/115 10 outliers
#in5=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.75))#78/115 20 outliers
#in6=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC"), nrep=20,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.75))#78/115 19 outliers
#in7=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=100,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#80/115 16 outliers &74
#plot_contaminated_fd(functional_nox_data(),ianos3$plot[2,],show_contam=T,split_class = T, split_contam = F)
#ino1=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#80/115 16 outliers &74
#ino2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))#80/115 12 outliers &68
#ino3=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.01,alphamin=c(0.85))#79/115 4 outliers &78
#ino4=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.85))#69/115 10 outliers &64
#ino6=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.5))#80/115 16 outliers &74
#ino7=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.75))#80/115 16 outliers &74
#ino5=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.9))#80/115 16 outliers &74
#ino8=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.5))#80/115 16 outliers &74
#ino9=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.75))#80/115 16 outliers &74
#ino10=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.9))#80/115 16 outliers &74
####in1=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=1,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#80/115 16 outliers &74
#in2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))#80/115 12 outliers &68
#in3=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.01,alphamin=c(0.85))#79/115 4 outliers &78
#in4=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.85))#69/115 10 outliers &64
#in6=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.5))#80/115 16 outliers &74
#in7=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.75))#80/115 16 outliers &74
#in5=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.9))#80/115 16 outliers &74
#in8=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.5))#80/115 16 outliers &74
#in9=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.75))#80/115 16 outliers &74
#in10=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC", "funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.9))#80/115 16 outliers &74


# ino1_fix=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.85)) #77/99
# ino12_fix=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))#78/99
# ino13_fix=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#71/99
# ino14_fix=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.01,alphamin=c(0.85))#63/99
 #ino5_fix=cluster_testdata(1, removal_nox(ino5), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.2,alphamin=c(0.85)) #80/103
#ino52_fix=cluster_testdata(1, removal_nox(ino5), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.85))#80/103
 #ino53_fix=cluster_testdata(1, removal_nox(ino5), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#63/103
 #ino54_fix=cluster_testdata(1, removal_nox(ino5), K=2, fun_names=c("funHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.01,alphamin=c(0.85))#77/103
#ino0_fix=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("cfunHDDC"), nrep=50,init="kmeans", simpre="NOx", threshold=0.1,alphamin=c(0.99999)) #78/99 8 outliers?
 removal_nox <- function(indata, nox = functional_nox_data()) {
  rem_ind <- which(indata$plot[2,]$contam == 0)
  nox$fd <- nox$fd[-rem_ind]
  nox$contam <- nox$contam[-rem_ind]
  nox$class <- nox$class[-rem_ind]
  return(nox)
}
plot_nox<-function(data1=functional_nox_data(nox, 8)){
  par(mfrow=c(1,3))
  plot(nox, col=CL[data1$class], main="a")
  rem_ind1 <- which(data1$class == 1)
  nox1=nox[-rem_ind1]
  plot(nox1, col=CL[2], main="b")
  rem_ind2 <- which(data1$class == 2)
  nox2=nox[-rem_ind2]
  plot(nox2, col=CL[1], main="c")
  
}
plot_nox_mic<-function(data1=functional_nox_data(nox, 8),c1=ino1$plot[2,]$class,
                       d1=t(ino1$plot[2,]$contam)){
  par(mfrow=c(1,3))
  plot(nox, col=CL[data1$class], main="a")
  c11<-ifelse(c1 == 1, 2, 1)
  func <- ifelse(d1 == 1, c11, 3)
  rem_ind1 <- which(c11 == 1)
  nox1=nox[-rem_ind1]
  func1=func[-rem_ind1]
  plot(nox1, col=CL[func1], main="b", ylab="")
  rem_ind2 <- which(c11 == 2)
  nox2=nox[-rem_ind2]
  func2=func[-rem_ind2]
  plot(nox2, col=CL[func2], main="c", ylab="")
}
plot_classnox<-function(data1=functional_nox_data(nox, 8), c1=ino1$plot[2,]$class,
                        c2=ino5$plot[2,]$class,d1=t(ino1$plot[2,]$contam),
                        d2=t(ino5$plot[2,]$contam)){
  c11<-ifelse(c1 == 1, 2, 1)
  func <- ifelse(d1 == 1, c11, 3)
  funco <- ifelse(d2 == 1, c2, 3)
  par(mfrow=c(2,2))
  rem_ind1 <- which(c11 == 1)
  nox1=nox[-rem_ind1]
  func1=func[-rem_ind1]
  plot(nox1, col=CL[func1], main="a",xlab="", ylab="")
  rem_ind2 <- which(c11 == 2)
  nox2=nox[-rem_ind2]
  func2=func[-rem_ind2]
  plot(nox2, col=CL[func2], main="b",xlab="", ylab="")
  rem_indo1 <- which(c2 == 1)
  nox1=nox[-rem_indo1]
  funco1=funco[-rem_indo1]
  plot(nox1, col=CL[funco1], main="c",xlab="", ylab="")
  rem_indo2 <- which(c2 == 2)
  nox2=nox[-rem_indo2]
  funco2=funco[-rem_indo2]
  plot(nox2, col=CL[funco2], main="d",xlab="", ylab="")
  
}
#inc1=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.5))#77/115 6 outliers
#inc2=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.75))#76/115 9 outliers
#inc3=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.85))#77/115 6 outliers
#inc4=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(0.9))#76/115 9 outliers
#inc5=cluster_testdata(1, removal_nox(ino1), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(1))#65/99 0 outliers
#inc6=cluster_testdata(1, removal_nox(ino5), K=2, fun_names=c("CNmixt"), nrep=50,init="kmeans", simpre="NOx", threshold=0.05,alphamin=c(1))#86/103 0 outliers
#in10=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC","funHDDC"), nrep=1,init="kmeans", simpre="NOx", threshold=0.45,alphamin=c(0.85))#94/115 etamx=100,all outliers &74
in11=cluster_testdata(1, functional_nox_data(), K=2, fun_names=c("cfunHDDC","funHDDC"), nrep=1,init="kmeans", simpre="NOx", threshold=0.8,alphamin=c(0.85))#94/115 etamx=100,all outliers &74

