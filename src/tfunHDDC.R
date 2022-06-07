library(stats)
library(funHDDC) # clustering library
# library(RColorBrewer) # color pallette library
#library(ContaminatedMixt)
# library(dplyr)
# library(tidyverse)
library(fda)
# library(rworldmap)
library(mclust)
library(foreach)
library(doFuture)
library(future.batchtools)

tfunHDDC  <-
  function(data, K=1:10, model="AkjBkQkDk",dfstart=50,dfupdate="approx",dconstr="no",threshold=0.1, itermax=200, eps=1e-6,init='random', criterion="bic",
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
      try(res <- .tfunhddc_main1(model=model, K=K, dfstart=dfstart,dfupdate=dfupdate,dconstr=dconstr,threshold=threshold, itermax=itermax, ...),silent=TRUE)
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
      if(n>1) cat("tfunHDDC: \n")
      
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

.tfunhddc_main1 <- function(fdobj, K, dfstart=dfstart,dfupdate=dfupdate,dconstr=dconstr,model, itermax=200,threshold, method, algo, eps=1e-6, init, init.vector, mini.nb, min.individuals, noise.ctrl, com_dim=NULL, kmeans.control, d_max, ...){
  
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
    tw<-matrix(0, N, K)
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
        prms <- .tfunhddc_main1(DATA, K,tfunhddc_main1, model, dfupdate,dconstr, threshold, method, algo, mini.nb[2], 0, 'random', mini.nb = mini.nb, min.individuals = min.individuals, noise.ctrl = noise.ctrl, com_dim = com_dim, d_max=d_max)
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
  } else {
    t <- matrix(1, N, 1) #@5 K IS 1
    tw<-matrix(1, N, 1)
  }
  
  nux <-c(rep.int(dfstart, K))
  initx<-.funhddt_init(fdobj, K,t,nux, model, threshold, method, noise.ctrl, com_dim, d_max)
  tw<-.funhddt_twinit(fdobj, initx,nux)
  
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
    
    m <- .funhddt_m_step1(fdobj, K, t,tw,nux,dfupdate,dconstr, model,threshold, method, noise.ctrl, com_dim, d_max) # A@5mstep is applied here
    #@A5 here we need to add vi and eta after t
    # print(m$mu)
    nux<-m$nux
    to <- .funhddt_e_step1(fdobj, m) #A@5 estep is applied here
    L <- to$L # A@5 this s the likelyhood  THIS HAS TO BE CHANGED
    t <- to$t#A@5 this is the new t
    tw<-to$tw
    #if (I>3)
    
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
    # print(paste0("lik=",likely[I] ))
    #print(paste0("test=",test ))
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
  nux <- matrix(m$nux, 1, m$K, dimnames=list(c(''), 'Degrees of freedom of t-distributions'=1:m$K))
  
  # A25 add alpha and eta
  
  # Other elements
  complexity <- .hdc_getComplexityt(m, p, dconstr)
  class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- 'hd' #A@5 alpha and eta will nedd class too
  cls <- max.col(t)#@ here the cluster is found
  
  ## we have to add the clasification as outlier or not
  converged = test<eps
  
  params = list(model=model, K=K, d=d, a=a, b=b, mu=mu, prop=prop,nux=nux, ev=m$ev, Q=m$Q,fpca=m$fpcaobj, loglik=likely[length(likely)], loglik_all = likely, posterior=t, class=cls, com_ev=com_ev, N=N, complexity=complexity, threshold=threshold, d_select=method, converged=converged) #A25 add alpha and eta
  #  params = list(model=model, K=K, d=d, a=a, b=b, mu=mu, prop=prop,alphax=alphax,etax=etax, ev=m$ev, Q=m$Q,fpca=m$fpcaobj, loglik=likely[length(likely)], loglik_all = likely, posterior=t, class=cls, com_ev=com_ev, N=N, complexity=complexity, threshold=threshold, d_select=method, converged=converged) #A25 add alpha and eta
  
  #A@5 CONTINUE FROM HERE
  
  
  # We compute the BIC / ICL
  bic_icl = .hdclassift_bic(params, p, dconstr)
  params$BIC = bic_icl$bic
  params$ICL = bic_icl$icl
  
  # We set the class
  class(params) <- 'tfunHDDC'
  
  return(params)
}
###############################################
#########################################################
.funhddt_init  <- function(fdobj, K,t,nux, model, threshold, method, noise.ctrl, com_dim, d_max){ #A@5 this is the init step
  if (class(fdobj)!='list') { x = t(fdobj$coefs) #A@5 THIS IS OUR CASE
  } else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
  #A@5 x is the coefficient in the fdobject
  # Some parameters
  N <- nrow(x)
  p <- ncol(x)
  prop <- c()
  n <- colSums(t) 
  prop <- n/N  #A5 props is pi as a vector needed for ai
  
  
  
  
  mu <- matrix(NA, K, p)
  
  # print(mu)
  ind <- apply(t>0, 2, which)
  n_bis <- c()
  for(i in 1:K) n_bis[i] <- length(ind[[i]])
  
  #
  #Calculation on Var/Covar matrices and mean mu vector
  #
  
  # we keep track of the trace (== sum of eigenvalues) to compute the b
  traceVect = c()
  
  
  ev <- matrix(0, K, p)
  Q <- vector(mode='list', length=K) #A@5 this is matrix q_k
  fpcaobj = list()
  for (i in 1:K){
    donnees <- .initmypca.fd1(fdobj,t[,i]) #A@5 we must change .mypca- bu everything else stays the same
    #A@% DON"T TOUCH
    #-----------------------------------
    mu[i,]<-donnees$mux
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
  
  
  
  list(model=model, K=K, d=d, a=ai, b=bi, mu=mu, prop=prop, ev=ev, Q=Q,fpcaobj = fpcaobj) #A@5 add the extra params(alpha and eta)
  
  
}

.initmypca.fd1 <- function(fdobj, Ti){
  if (class(fdobj)=="list"){ #A@5 NOT OUR CASE 
    #sauvegarde de la moyenne avant centrage
    mean_fd<-list()
    for (i in 1:length(fdobj)){
      mean_fd[[i]]<-fdobj[[i]]
    }
    
    #Cr?ation de la matrice de coef
    coef<-t(fdobj[[1]]$coefs)
    for (i in 2:length(fdobj)){
      coef<-cbind(coef,t(fdobj[[i]]$coefs))
    }
    coef1<-t(fdobj[[1]]$coefs)
    wtcov <- cov.wt(coef1,wt=Ti,method="ML")
    mat_cov <-wtcov$cov
    coefmean <- wtcov$center
    fdobj[[1]]$coefs <- sweep(fdobj[[1]]$coefs, 1, coefmean)  
    mean_fd[[1]]$coefs = as.matrix(data.frame(mean=coefmean))
    for (i in 2:length(fdobj)){ 
      coef1<-t(fdobj[[i]]$coefs)
      wtcov <- cov.wt(coef1,wt=Ti,method="ML")
      mat_cov <-rbind(mat_cov,wtcov$cov)
      coefmean <- wtcov$center
      fdobj[[i]]$coefs <- sweep(fdobj[[i]]$coefs, 1, coefmean)  
      mean_fd[[i]]$coefs = as.matrix(data.frame(mean=coefmean))
    }
    #Matrice de covariance #A@5 mat_cov is formula 19 without the big bracket which we have to add
    
    #Centrer les objets fonctionnels par groupe
    
    
    
    
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
    
    
    
    
    #mat_interm<-1/sqrt((sum(Ti)-1))*coef%*%(W_tot)^(1/2)
    
    #D?but Partie mise en commentaire
    #mat_interm<-1/sqrt((sum(Ti)))*coef%*%chol(W_tot,pivot=TRUE)
    #cov<-(.repmat(Ti,n=pcol,p=1)*t(mat_interm))%*%mat_interm
    #Fin de partie mise en commentaire
    
    #D?but correction
    #Construction matrice triangulaire de Choleski
    W_m <-  chol(W_tot) #A@5 this is sqrt(W)
    #Matrice de covariance 
    
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
    coef<-t(fdobj$coefs)
    wtcov <- cov.wt(coef,wt=Ti,method="ML")
    coefmean <- wtcov$center
    
    #Matrice de covariance #A@5 mat_cov is formula 19 without the big bracket which we have to add
    mat_cov <-wtcov$cov
    #Centrer les objets fonctionnels par groupe
    
    
    fdobj$coefs <- sweep(fdobj$coefs, 1, coefmean) #A@5 what does sweep do?substracts coefmean from gamma, i.e gammai-mui
    mean_fd$coefs = as.matrix(data.frame(mean=coefmean))
    
    #Calcul de la matrice des produits scalaires #A@5 this stays the SAME CAN"T TOUCH THIS!
    #----------------------------------------------------------
    W<-inprod(fdobj$basis,fdobj$basis)
    #pour ?viter les soucis num?riques on arrondit ? 0 les produits scalaires tr?s petits
    W[W<1e-15]=0
    #on centre les coefficients
    
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
    
    ipcafd <-list(valeurs_propres=valeurs_propres,harmonic=fonctionspropres,scores=scores,covariance=cov,U=bj,meanfd=mean_fd,mux=coefmean, W=W)
    #-------------------------------------------    
  }
  class(ipcafd) <- "ipca.fd"
  return(ipcafd)
}
###################################################
############################################################
.funhddt_twinit  <- function(fdobj, par,nux){ #A@5 this is the E step
  if (class(fdobj)!='list') {x = t(fdobj$coefs)}  #A@5 THIS IS OUR CASE
  else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))} #NOT THIS
  p <- ncol(x)
  N <- nrow(x)
  K <- par$K
  a <- par$a
  b <- par$b
  mu <- par$mu
  d <- par$d
  Q <- par$Q
  W <- matrix(0,N,K)
  
  
  
  #A@5 add the alhpas and the etas
  
  b[b<1e-6] <- 1e-6
  
  #############################aici am ajuns
  #################################################
  
  mah_pen<-matrix(0,N,K)
  
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
    mah_pen[,i]=.imahalanobis(x, muki,Wki ,Qk, aki, bki)
    W[,i] <- (nux[i]+p)/(nux[i]+mah_pen[,i])
  }
  return(W)
}
######################################################################
############################################################
.funhddt_e_step1  <- function(fdobj, par){ #A@5 this is the E step
  if (class(fdobj)!='list') {x = t(fdobj$coefs)}  #A@5 THIS IS OUR CASE
  else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))} #NOT THIS
  p <- ncol(x)
  N <- nrow(x)
  K <- par$K
  nux<-par$nux
  a <- par$a
  b <- par$b
  mu <- par$mu
  d <- par$d
  prop <- par$prop
  Q <- par$Q
  
  #A@5 add the alhpas and the etas
  
  b[b<1e-6] <- 1e-6
  ##################################################
  #########################################################
  ###########################################
  t <- matrix(0,N,K)
  tw<-matrix(0,N,K)
  mah_pen<-matrix(0,N,K)
  K_pen <- matrix(0,N,K)
  num <- matrix(0,N,K)
  ft <- matrix(0,N,K)
  #############################aici am ajuns
  #################################################
  
  
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
    mah_pen[,i]=.imahalanobis(x, muki,Wki ,Qk, aki, bki)
    tw[,i] <- (nux[i]+p)/(nux[i]+mah_pen[,i])
    K_pen[,i]<-log(prop[i])+lgamma((nux[i]+p)/2)-(1/2)*(s[i]+(p-d[i])*log(b[i]))-
      ((p/2)*(log(pi)+log(nux[i]))+lgamma(nux[i]/2)+((nux[i]+p)/2)*
         (log(1+mah_pen[,i]/nux[i])))
    
  }
  ft <- exp(K_pen)
  #ft_den used for likelihood
  ft_den <- rowSums(ft)
  kcon <- -apply(K_pen,1,max)
  K_pen <- K_pen + kcon
  num <- exp(K_pen)
  t <- num/rowSums(num)
  ###############################################################tw#######
  
  
  
  
  #pdf_pen=t(prop*Stot_pen)
  #############################################
  
  
  
  # vx[which(Stot_pen<1e-10)]=0
  
  #what is this L?
  L1<- sum(log(ft_den))
  L<-sum(log(rowSums(exp(K_pen)))-kcon)
  # print(paste0("log check ",L1-L ))
  # A <- -1/2*t(K_pen) 
  #L<-sum(log(rowSums(ax_pen*exp(A-apply(A,1,max))))+apply(A,1,max))
  #likelihood as in CNmixt
  # L<-sum(log(rowSums(pdf_pen)))
  #A25 add calculations for v (formula 13)exp
  
  # t <- matrix(0,N,K)
  # for (i in 1:K) t[,i] <- ax_pen[,i]/rowSums(t(exp(-0.5*K_pen+0.5*matrix(K_pen[i,],K,N,byrow = TRUE)))*ax_pen) #A@5 this  1 over is the exponential part on the denominator on formula 12  (the inverse)
  #for (i in 1:K) t[,i] <- ax_pen[,i]*exp(A[,i]-apply(A,1,max))/rowSums(ax_pen*exp(A-apply(A,1,max)))
  
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
  list(t=t,tw=tw,L=L) #A@5 add v to t and return it  as well
}

# A@5 pass the z, v, and eta in a names list from the e step par

.funhddt_m_step1  <- function(fdobj, K, t,tw,nux,dfupdate, dconstr,model, threshold, method, noise.ctrl, com_dim, d_max){ #A@5 this is the M step
  if (class(fdobj)!='list') { x = t(fdobj$coefs) #A@5 THIS IS OUR CASE
  } else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
  #A@5 x is the coefficient in the fdobject
  # Some parameters
  N <- nrow(x)
  p <- ncol(x)
  prop <- c()
  n <- colSums(t) 
  prop <- n/N  #A5 props is pi as a vector
  
  #print(paste0("calcalphax=",alphax ))
  
  
  
  
  mu <- matrix(NA, K, p)
  mu1 <- matrix(NA, K, p)
  corX=matrix(0,N,K)
  corX=t*tw
  for (i in 1:K) {
    #mu[i, ] <- colSums(x*t[, i])/n[i] #A@5 add the big parenthesis from formula 18 and becomes something like colSums(x*t[,i]* a slive of v)
    #A@5 and n[i] must be multiplied by some form of colsum of the product between t and the same thing as above
    # a@5 new
    
    mu[i,]<-apply(t(as.matrix(corX[,i]) %*% matrix(1,1,p)) * t(x), 1, sum) / sum(corX[,i])
    mu1[i,]<-colSums(corX[,i]*x)/sum(corX[,i])
    
    #mu[i,]<-colSums(x*matrix(corX[,i],N,p,byrow=FALSE))/sum(corX[,i])
  }
  #  print(max(mu-mu1))
  # print(mu)
  # for (i in 1:K) mu[i, ] <- colSums((z[, i]*(v[, i] + (1-v[, i])/eta[i]))*x)/colSums((z[, i]*(v[, i] + (1-v[, i])/eta[i])))
  ind <- apply(t>0, 2, which)
  n_bis <- c()
  for(i in 1:K) n_bis[i] <- length(ind[[i]])
  #calculation of the degree
  if(dfupdate=="approx"){
    testing <- try(jk861 <- tyxf8(dconstr,nux,n,t,tw,K,p,N),silent=TRUE)
    if(all(is.finite(testing))){
      nux <- jk861
    }
    else{break}
  }
  else{
    if(dfupdate=="numeric"){
      testing <- try(jk861 <- tyxf7(dconstr,nux,n,t,tw,K,p,N),silent=TRUE)
      if(all(is.finite(testing))){
        nux <- jk861
      }
      else{break}
    }
  }
  
  
  #
  #Calculation on Var/Covar matrices
  #
  
  # we keep track of the trace (== sum of eigenvalues) to compute the b
  traceVect = c()
  
  
  ev <- matrix(0, K, p)
  Q <- vector(mode='list', length=K) #A@5 this is matrix q_k
  fpcaobj = list()
  for (i in 1:K){
    donnees <- .mypcat.fd1(fdobj,t[,i],corX[,i]) #A@5 same as for CfunHDDC
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
  
  
  
  #A@5 add the CM2 step
  
  # A@5 new mahalanobis rowSums(A^2)+1/b[i]*rowSums(B^2) from e step
  ###########new mu in terms of new eta ca in CNMixt
  
  ###################################
  #A@5 check the other library on how to use optimize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  list(model=model, K=K, d=d, a=ai, b=bi, mu=mu, prop=prop, nux=nux,ev=ev, Q=Q,fpcaobj = fpcaobj) #A@5 add the extra params(alpha and eta)
  
  
}
#############################################
##################################################
#degree functions
tyxf7 <-
  function(dconstr,nux,n,t,tw,K,p,N){
    if(dconstr=="no"){
      dfoldg <- nux
      for(i in 1:K){
        constn <- 1 + (1/n[i]) * sum(t[,i]*(log(tw[,i])-tw[,i])) + digamma((dfoldg[i]+p)/2)-log((dfoldg[i]+p)/2)
        ## 			dhfgs78[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn ,
        ## 					lower=0.0001, upper=1000, tol=0.00001)$root
        
        nux[i] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn ,
                           lower=0.0001, upper=1000, tol=0.00001)$root
        ##constn <- -constn
        ##print(constn)
        ##print(paste("fix", dhfgs78[g] - ((-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn)))))
        ##print(paste("old", dhfgs78[g] - (-exp(constn)/(1-exp(constn)))))
        ##print(c((2-sqrt(4+4*4*exp(constn)))/2, (2+sqrt(4+4*4*exp(constn)))/2))
        ##print((-exp(constn))/(1-exp(constn)))
        if(nux[i]>200){
          nux[i]<-200
        }
        if(nux[i]<2){
          nux[i]<-2
        }
      }
    }
    else
    {
      dfoldg <- nux[1]
      constn <- 1 + (1/N) * sum(t * (log(tw)-tw)) + digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
      
      
      dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn,
                             lower=0.0001, upper=1000, tol=0.01)$root
      if(dfsamenewg>200){
        dfsamenewg<-200
      }
      if(dfsamenewg<2){
        dfsamenewg<-2
      }
      nux <- c(rep(dfsamenewg,K))
    }
    return(nux)
  }

tyxf8 <- function(dconstr,nux,n,t,tw,K,p,N){
  if(dconstr=="no"){
    dfoldg <- nux
    for(i in 1:K){
      constn <- 1 + (1/n[i]) * sum(t[,i]*(log(tw[,i])-tw[,i])) + digamma((dfoldg[i]+p)/2) 
      - log((dfoldg[i]+p)/2)
      ## 			dhfgs78[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn ,
      ## 					lower=0.0001, upper=1000, tol=0.00001)$root
      constn <- -constn
      ##print(constn)
      ##print(paste("fix", dhfgs78[g] - ((-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn)))))
      ##print(paste("old", dhfgs78[g] - (-exp(constn)/(1-exp(constn)))))
      ##print(c((2-sqrt(4+4*4*exp(constn)))/2, (2+sqrt(4+4*4*exp(constn)))/2))
      ##print((-exp(constn))/(1-exp(constn)))
      nux[i] <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[i]/2))-(dfoldg[i]/2-1/2)))/(1-exp(constn))
      ##     dhfgs78[g] <- (-exp(constn))/(1-exp(constn))
      if(nux[i]>200){
        nux[i]<-200
      }
      if(nux[i]<2){
        nux[i]<-2
      }
    }
  }
  else
  {
    dfoldg <- nux[1]
    constn <- 1 + (1/N) * sum(t * (log(tw)-tw)) + digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
    constn <- -constn
    ##		dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn,
    ##				lower=0.0001, upper=1000, tol=0.01)$root
    dfsamenewg <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg/2))-(dfoldg/2-1/2)))/(1-exp(constn))
    if(dfsamenewg>200){
      dfsamenewg<-200
    }
    if(dfsamenewg<2){
      dfsamenewg<-2
    }
    nux <- c(rep(dfsamenewg,K))
  }
  return(nux)
}
#######################################################
############################################################
.mypcat.fd1 <- function(fdobj, Ti,CorI){
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

.hdclassift_bic <- function(par, p, dconstr,data=NULL){
  model <- par$model
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  N <- par$N
  prop <- par$prop
  mux<-par$mux
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
    mah_pen<-matrix(0, N, K)
    K_pen <- matrix(0, K, N)
    ft<-matrix(0, K, N)
    for (i in 1:K) {
      s <- sum(log(a[i, 1:d[i]]))
      X <- data-matrix(mu[i, ], N, p, byrow=TRUE)
      proj <- (X%*%Q[[i]])%*%t(Q[[i]])
      A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i, 1:d[i]], d[i]))
      B <- X-proj
      mah_pen[,i]<-rowSums(A^2)+1/b[i]*rowSums(B^2)
      K_pen[,i]<-log(prop[i])+lgamma((nux[i]+p)/2)-(1/2)*(s+(p-d[i])*log(b[i]))-
        ((p/2)*(log(pi)+log(nux[i]))+lgamma(nux[i]/2)+((nux[i]+p)/2)*
           (log(1+mah_pen[,i]/nux[i])))
    }
    ft <- exp(K_pen)
    #ft_den used for likelihood- to erese one of them L or L1
    ft_den <- rowSums(ft)
    kcon <- -apply(K_pen,1,max)
    K_pen <- K_pen + kcon
    L<-sum(log(rowSums(exp(K_pen)))-kcon)
    L1<- sum(log(ft_den))
    ###############################################################tw#######
    
    #what is this L?
  } else L <- par$loglik[length(par$loglik)]
  
  #@ add K or 1 parameters for nux in ro
  if (dconstr=="no")
    ro <- K*(p+1)+K-1
  else
    ro <- K*p+K
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

.hdc_getComplexityt = function(par, p, dconstr){
  model <- par$model
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  
  prop <- par$prop
  #@ add in ro K parameters for nuk or add 1 parameter if degrees equal
  if (dconstr=="no")
    ro <- K*(p+1)+K-1
  else
    ro <- K*p+K
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
