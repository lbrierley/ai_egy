#### EMBARCADERO FUNCTIONS ADJUSTED TO ALLOW CUSTOM BART PARAMETERS: ####

bart.flex <- function(x.data, y.data,
                      x.weight = NULL,
                      y.name = NULL,
                      n.trees = 200,
                      k = 2.0, power = 2.0, base = 0.95) {
  
  train <- cbind(y.data, x.data) 
  if(!is.null(y.name)) {colnames(train)[1] <- y.name}
  train <- na.omit(train)
  model <- bart(y.train = train[,1], 
                x.train = train[,2:ncol(train)], 
                weights = x.weight,
                k = k,
                power = power,
                base = base,
                ntree = n.trees,
                keeptrees=TRUE,
                nchain = N_CHAINS,
                nthread = N_CHAINS)
  return(model)
}

varimp.diag <- function(x.data, y.data, x.weight = NULL, iter=50, quiet=FALSE,
                        k = 2.0, power = 2.0, base = 0.95) {
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist <- colnames(x.data)
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY 4 THIS CODE :) 
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data,
                               x.weight = x.weight,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  for (n.trees in c(10, 20, 50, 100, 150, 200)) {
    
    cat(paste('\n', n.trees, 'tree models:', iter, 'iterations\n'))
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   x.weight = x.weight,
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      vi.j <- varimp(model.j)
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    
    if(n.trees==10) { vi <- vi.j } else {  vi <- cbind(vi,vi.j[,2])  }
  }
  
  colnames(vi) <- c('variable','10','20','50','100','150','200')
  vi <- reshape::melt(vi, "variable")
  colnames(vi) <- c('variable','trees','imp')
  
  vi %>% group_by(variable) %>% summarise(max = max(imp)) -> vi.fac
  vi.fac <- vi.fac[order(-vi.fac$max),]
  
  vi$names <- factor(vi$variable, levels=vi.fac$variable)
  
  g1 <- ggplot2::ggplot(vi, aes(y=imp, x=names, group=trees)) +
    geom_line(aes(color=trees)) + geom_point(size=3) + theme_classic() +
    ylab("Relative contribution\n") + xlab("\nVariables dropped") +
    ggpubr::rotate_x_text(angle = 35) + 
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=14,face="bold")); print(g1)
  
}

variable.step <- function(x.data, y.data, x.weight = NULL, n.trees=10, iter=50, quiet=FALSE,
                          k = 2.0, power = 2.0, base = 0.95) {
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  comp <- complete.cases(x.data)
  
  if(length(comp) < (nrow(x.data))) {
    message("Some rows with NA's have been automatically dropped. \n")
  }
  x.data <- x.data[comp,]
  y.data <- y.data[comp]
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               x.weight = x.weight,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames) > 0) {
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-any_of(dropnames)) -> x.data  
  
  ###############
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist.orig <- varlist <- colnames(x.data)
  
  rmses <- data.frame(Variable.number=c(),RMSE=c())
  dropped.varlist <- c()
  
  for(var.j in c(nvars:3)) {
    
    print(noquote(paste("Number of variables included:",var.j)))
    print(noquote("Dropped:"))
    print(if(length(dropped.varlist)==0) {noquote("")} else {noquote(dropped.varlist)})
    
    rmse.list <- c()
    
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   x.weight = x.weight,
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      quietly(vi.j <- varimp(model.j))
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      
      pred.p <- colMeans(pnorm(model.j$yhat.train))[y.data==1]
      pred.a <- colMeans(pnorm(model.j$yhat.train))[y.data==0]
      
      pred.c <- c(pred.p, pred.a)
      true.c <- c(rep(1,length(pred.p)), rep(0,length(pred.a)))
      rmsej.i <- Metrics::rmse(true.c,pred.c)
      rmse.list <- c(rmse.list,rmsej.i)
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    vi.j <- vi.j[order(vi.j[,2]),]
    
    drop.var <- vi.j[1,1]
    dropped.varlist <- c(dropped.varlist,as.character(drop.var))
    
    rmsej <- mean(rmse.list)
    
    rmses <- rbind(rmses,c(nvars-var.j,rmsej)); colnames(rmses) <- c('VarsDropped','RMSE')
    
    varnums <- varnums[!(varnums==which(varlist.orig==drop.var))]
    varlist <- varlist.orig[varnums]
    print(noquote("---------------------------------------"))
  }
  
  g1 <- ggplot2::ggplot(rmses, aes(y=RMSE, x=VarsDropped)) +
    geom_line(color="black") + geom_point(size=3) + theme_bw() +
    ylab("RMSE of model\n") + xlab("\nVariables dropped") +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14,face="bold")) +
    scale_x_discrete(limits=c(0:(nrow(rmses)))); print(g1)
  
  print(noquote("---------------------------------------"))
  print(noquote("Final recommended variable list"))
  varlist.final <- varlist.orig[!(varlist.orig %in% dropped.varlist[0:(which(rmses$RMSE==min(rmses$RMSE))-1)])]
  print(noquote(varlist.final))
  invisible(varlist.final)
}

bart.step <- function(x.data, y.data, x.weight = NULL,
                      iter.step=100, tree.step=10,
                      iter.plot=100,
                      k = 2.0, power = 2.0, base = 0.95,
                      full=FALSE,
                      quiet=FALSE) {
  
  ###############
  
  # auto-drops 
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               x.weight = x.weight,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  quiet2 <- quiet
  if(full==TRUE){varimp.diag(x.data, y.data, iter=iter.plot, x.weight, quiet=quiet2)}
  vs <- variable.step(x.data,
                      y.data,
                      x.weight=x.weight,
                      n.trees=tree.step,
                      iter=iter.step,
                      quiet=quiet2)
  
  invisible(best.model <- bart.flex(x.data = x.data[,vs], y.data = y.data, 
                                    x.weight = x.weight,
                                    n.trees=200,
                                    k = k,
                                    power = power,
                                    base = base))
  if(full==TRUE){varimp(best.model, plots=TRUE)}
  if(full==TRUE) {p <- summary(best.model, plots=TRUE)
  print(p)} else 
  {p <- summary(best.model, plots=FALSE)
  print(p)}
  invisible(best.model)
}

#### Next two functions are adaptations of stuff for calculating metrics done inside embarcadero::summary ####

get_threshold <- function(object){
  fitobj <- object$fit
  
  true.vector <- fitobj$data@y 
  
  pred <- prediction(colMeans(pnorm(object$yhat.train)), true.vector)
  
  perf.tss <- performance(pred,"sens","spec")
  tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
  tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
  return(thresh)
}

get_sens_and_spec <- function(sdm, xtest, ytest, cutoff){
  pred <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm, type = "bart") %>%
    pnorm() %>%
    colMeans() %>%
    prediction(labels = ytest[which(complete.cases(xtest))])
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  sens <- perf@x.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  spec <- perf@y.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  tss <- tss_df[which.min(abs(tss_df$alpha-cutoff)),'tss']
  auc <- performance(pred,"auc")@y.values[[1]]
  return(pairlist("pred"=pred,
                  "perf"=perf,
                  "sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "auc"=auc))
}

get_sens_and_spec_for_single_model <- function(pred){
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  # cutoff <- min(tss_df$alpha[which(tss_df$tss==max(tss_df$tss))])
  sens <- perf@x.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  spec <- perf@y.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  tss <- tss_df[which.min(abs(tss_df$alpha-cutoff)),'tss']
  auc <- performance(pred,"auc")@y.values[[1]]
  return(pairlist("sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "auc"=auc))
}

get_sens_and_spec_ci <- function(sdm, xtest, ytest, ri, cutoff){
  predmat <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm, type = "bart", group.by=ri) %>%
    pnorm()
  pred_by_model <- sapply(1:nrow(predmat),
                          FUN=function(i){predmat[i,] %>%
                              prediction(labels = ytest[which(complete.cases(xtest))])})
  metrics_by_model <- lapply(1:length(pred_by_model),
                             FUN=function(i){
                               get_sens_and_spec_for_single_model(pred_by_model[[i]])})
  sens_by_model <- sapply(1:length(metrics_by_model),
                          FUN=function(i){metrics_by_model[[i]]$sens})
  spec_by_model <- sapply(1:length(metrics_by_model),
                          FUN=function(i){metrics_by_model[[i]]$spec})
  tss_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$tss})
  auc_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$auc})
  return(pairlist("sens"=quantile(sens_by_model, c(.025, .975), names=FALSE),
                  "spec"=quantile(spec_by_model, c(.025, .975), names=FALSE),
                  "tss"=quantile(tss_by_model, c(.025, .975), names=FALSE),
                  "auc"=quantile(auc_by_model, c(.025, .975), names=FALSE)))
}

