############################
#### computation of MDL ####
############################

# chrom0stat is a list:
# chrom0stat$chrom
# chrom0stat$loglike
# chrom0stat$mdl.pen

MDL <- function(chrom, x1, x2, y, A, delta.t, delta.w, opt=get.opt.spec(),
                min.span=5, assign.emiss=NA, emdl=TRUE, v1=1, v2=1, chrom0stat=NULL){
  minx1 <- min(x1)-1
  maxx1 <- max(x1)+1

  B <- sum(chrom)+1
  x1.index <- cut(x1,breaks=c(minx1,x1[c(FALSE,chrom)],maxx1),labels=(1:B),right=FALSE)
  x1.len <- tapply(x1,x1.index,length)
  if (min(x1.len)>=min.span){
    mdl.pen <- array(dim=B)
    loglike <- array(dim=B)

    ## use chrom0stat or not
    if (is.null(chrom0stat)){
      ks <- 1:B
    } else {
      cres <- compare.chroms(chrom, chrom0stat$chrom)
      ks <- which(is.na(cres))
      for (j in which(!is.na(cres))){
        mdl.pen[j] <- chrom0stat$mdl.pen[cres[j]]
        loglike[j] <- chrom0stat$loglike[cres[j]]
      }
    }

    for (k in ks){
      tempy <- apply(y[,x1.index==k],1,sum)
      tempfit <- spec(x=x2, y=tempy, A=A, delta.t=delta.t, delta.w=delta.w,reps=x1.len[k],display=F,opt=opt,assign.emiss=assign.emiss,emdl=emdl,v=v1,simple=T)
      if (is.infinite(tempfit$mdl.pen[tempfit$mdl.ind[1],tempfit$mdl.ind[2]])){
        return(list(mdl=Inf))
      }
      loglike[k] <- sum(dpois(as.vector(y[,x1.index==k]), lambda=rep(exp(tempfit$off.const+tempfit$fitted.logbright-log(x1.len[k])),x1.len[k]),log=T))
      mdl.pen[k] <- tempfit$mdl.pen[tempfit$mdl.ind[1],tempfit$mdl.ind[2]]
    }
    #-sum(loglike)+sum(mdl.pen)++B*log(n1)+sum(log(x1.len))
    list(mdl=-sum(loglike)+sum(mdl.pen)*v1+(log(B)+sum(log(x1.len)))*v2,
         loglike=loglike, mdl.pen=mdl.pen)
  } else {
    list(mdl=Inf)
  }
}


# detect the same pieces in chrom0
compare.chroms <- function(chrom, chrom0){
  bs <- c(1, which(chrom), length(chrom)+1)
  bs0 <- c(1, which(chrom0), length(chrom0)+1)
  B <- sum(chrom)+1
  matchid <- array(dim=B)
  for (i in (1:B)){
    j <- which(bs0==bs[i])
    if (length(j)>0){
      if (bs[i+1] == bs0[j+1]){
        matchid[i] <- j
      }
    }
  }
  return(matchid)
}

###################################################################
#### growing trees algorithm for detecting changes across time ####
###################################################################
# opt is previously set as get.opt.spec(nlambda=200,lam.ratio.len=200)

spec.tbreak <- function(x1, x2, y, A, delta.t, delta.w,
                        display=TRUE, opt=get.opt.spec(),
                        cpus=4, assign.emiss=NA, emdl=TRUE, v1=1, v2=1,
                        simple=TRUE, max.B=4, min.span=5){
  if (cpus>1){
    cl <- parallel::makeCluster(min(cpus, parallel::detectCores()))
    doParallel::registerDoParallel(cl)
  }

  n1 <- length(x1)
  objlist <- list(list())

  nbreaks <- 1
  cont <- T
  breaks.ind <- ((min.span+1):(n1-min.span+1))-1
  cbreaks <- c()
  #old.min.criteria <- Inf
  res <- MDL(chrom=rep(F,(n1-1)),x1=x1,x2=x2,y=y,A=A,delta.t=delta.t,delta.w=delta.w,opt=opt,min.span=min.span,assign.emiss=assign.emiss,emdl=emdl, v1=v1, v2=v2, chrom0stat=NULL)
  old.min.criteria <- res$mdl
  best.criteria.seq <- old.min.criteria
  chrom0stat <- list(chrom=rep(F,(n1-1)), loglike=res$loglike, mdl.pen=res$mdl.pen)
  while ((nbreaks<=(max.B-1))&&cont){

    S <- length(breaks.ind)
    pop <- matrix(F,nrow=S,ncol=(n1-1))
    pop[,cbreaks] <- T
    for (i in (1:S)){
      pop[i,breaks.ind[i]] <- T
    }

    cat(nbreaks,"\n")
    # compute mdl
    if (cpus>1){
      #criteria
      res <- foreach::foreach(mci = 1:S,.verbose=display) %dopar% {
        MDL(chrom=pop[mci,],x1=x1,x2=x2,y=y,A=A,delta.t=delta.t,delta.w=delta.w,opt=opt,min.span=min.span,assign.emiss=assign.emiss,emdl=emdl, v1=v1, v2=v2, chrom0stat=chrom0stat)
      }
    } else {
      #criteria <- array(dim=S)
      res <- list()
      for (mci in (1:S)){
        #criteria[mci]
        res[[mci]] <- MDL(chrom=pop[mci,],x1=x1,x2=x2,y=y,A=A,delta.t=delta.t,delta.w=delta.w,opt=opt,min.span=min.span,assign.emiss=assign.emiss,emdl=emdl, v1=v1, v2=v2, chrom0stat=chrom0stat)
      }
    }

    criteria <- sapply(res, function(x){x$mdl})
    objlist[[nbreaks]] <- list(nbreaks=nbreaks,pop=pop,criteria=criteria)

    # updating
    min.ind <- which.min(criteria)
    min.criteria <- min(criteria)
    best.criteria.seq <- c(best.criteria.seq, min.criteria)
    if (min.criteria < old.min.criteria){
      #best.obj <- objlist[[nbreaks]]
      cbreaks <- c(cbreaks,breaks.ind[min.ind])
      nbreaks <- nbreaks + 1
      breaks.ind <- ((min.span+1):(n1-min.span+1))-1
      for (i in (1:length(cbreaks))){
        breaks.ind <- setdiff(breaks.ind,((cbreaks[i]-min.span+1):(cbreaks[i]+min.span-1)))
      }
      if (length(breaks.ind)==0){
        cont <- F
      }
      old.min.criteria <- min.criteria
    } else {
      cont <- F
    }

    # construct chrom0stat for next iteration
    if (cont){
      chrom0 <- (rep(F,n1-1)); chrom0[cbreaks] <- TRUE
      chrom0stat <- list(chrom=chrom0, loglike=res[[min.ind]]$loglike, mdl.pen=res[[min.ind]]$mdl.pen)
    }
  }
  if (cpus>1) parallel::stopCluster(cl)

  # final fit
  minx1 <- min(x1)-1
  maxx1 <- max(x1)+1
  best.chrom <- (rep(FALSE,n1-1))
  best.chrom[cbreaks] <- TRUE
  B <- sum(best.chrom)+1
  x1.index <- cut(x1,breaks=c(minx1,x1[c(FALSE,best.chrom)],maxx1),labels=(1:B),right=FALSE)
  x1.len <- tapply(x1,x1.index,length)
  best.fit <- array(list(),dim=B)
  for (k in (1:B)){
    tempy <- apply(y[,x1.index==k],1,sum)
    best.fit[[k]] <- spec(x=x2, y=tempy, A=A, delta.t=delta.t, delta.w=delta.w,reps=x1.len[k],display=F,opt=opt,assign.emiss=assign.emiss,emdl=emdl,v=v1,simple=simple)
  }
  return(list(B=B,best.fit=best.fit,breaks=c(1,cbreaks+1),content=objlist,
              best.criteria.seq=best.criteria.seq, opt=opt,
              assign.emiss=assign.emiss, emdl=emdl, v1=v1, v2=v2, max.B=max.B,
              min.span=min.span))
}


########################
#### heatmap matrix ####
########################
# fit corresponds to out$best.fit
# breaks corresponds to out$breaks

create.specheatmap.matrix <- function(x1, x2, breaks, fit, log=FALSE){
  B <- length(breaks)
  temp <- matrix(nr=length(x2), nc=length(x1))
  tbreaks <- sort(c(breaks,length(x1)+1))
  for (j in (1:B)){
    temp[,tbreaks[j]:(tbreaks[j+1]-1)] <- (fit[[j]]$fitted.logbright)
  }
  dimnames(temp) <- list(x2, x1)
  if (log){
    return(t(temp))
  } else {
    return(t(exp(temp)))
  }
}
