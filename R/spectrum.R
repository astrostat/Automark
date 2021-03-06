##################################
#### default options for spec ####
##################################
get.opt.spec <- function(p=3, K=30, alpha=1, nlambda=100, bs="rb", lam.ratio.lims=c(1e-5,1-1e-5),
                         lam.ratio.len=100, glmnet.lam.min=0.001, bgamma=1){
  opt <- list(p=p,
              K=K,
              alpha=alpha,
              nlambda=nlambda,
              bs=bs,
              lam.ratio.lims=lam.ratio.lims,
              lam.ratio.len=lam.ratio.len,
              glmnet.lam.min=glmnet.lam.min,
              bgamma=bgamma)
  return(opt)
}

##########################
#### spectrum fitting ####
##########################

spec <- function(x, y, A, delta.t, delta.w, reps=1, display=TRUE,
                 opt=get.opt.spec(), assign.emiss=NA, emdl=TRUE, simple=TRUE, v=1){

  # for assign.emiss=NULL, nlambda should be greatly increased
  # it will be automatically set to 10000 if it is set smaller
  # if this change is not desired, one can set lam.ratio.len=-1
  # note: lam.ratio.len is useless in the setting of assign.emiss=NULL
  if ((is.null(assign.emiss))&&(opt$lam.ratio.len!=-1)){
    if (opt$nlambda<10000) opt$nlambda <- 10000
  }

  # assign options
  p<-opt$p
  K<-opt$K
  alpha<-opt$alpha
  nlambda <-opt$nlambda
  bs<-opt$bs
  lam.ratio.lims<-opt$lam.ratio.lims
  lam.ratio.len<-opt$lam.ratio.len
  glmnet.lam.min<-opt$glmnet.lam.min
  bgamma<-opt$bgamma


  # main
  n <- length(y)
  off.const <- log(delta.t)+log(delta.w)+log(A)+log(reps)
  if (!is.null(assign.emiss)){
    if (is.na(assign.emiss)) assign.emiss <- 1:n
    assign.emiss <- sort(assign.emiss)
  }


  # select splines
  if (bs=="tr"){
    basis <- smooth.construct.tr.smooth.spec(mgcv::s(x, bs="tr", m=p, k=p+K+1) ,data.frame(cbind(x,y)),NULL)
    B.splines <- basis$X[,-1]
  } else if (bs=="rb"){
    basis <- smooth.construct.rb.smooth.spec(mgcv::s(x, bs="rb", m=p, k=p+K+1) ,data.frame(cbind(x,y)),NULL)
    B.splines <- basis$X[,-1]
  } else if (bs=="bs"){
    knot<-quantile(x[2:(n-1)],seq(0,1,length=K+2))[2:(K+1)]
    B.splines <- as.matrix(splines::bs(x,knots=knot,degree=p,intercept=F))
  }

  # create data matrix for atoms
  n1 <- ncol(B.splines)
  if (is.null(assign.emiss)){
    n2 <- 0
    X <- B.splines
    lam.ratios <- 1
    lam.ratio.len <- 1
  } else {
    B.atoms <- as.matrix(diag(rep(1,n))[,assign.emiss])
    n2 <- ncol(B.atoms)
    X <- cbind(B.splines,B.atoms)
    lam.ratios <- seq(lam.ratio.lims[1],lam.ratio.lims[2],len=lam.ratio.len)
  }

  # for ratios
  #lam.ratios <- exp(seq(log(lam.ratio.lims[1]),log(lam.ratio.lims[2]),len=lam.ratio.len))
  mdl <- matrix(nrow=lam.ratio.len,ncol=nlambda)
  mdl.pen <- matrix(nrow=lam.ratio.len,ncol=nlambda)
  beta.cube <- array(dim=c(lam.ratio.len,nlambda,n1+n2+1))
  lams <- matrix(nrow=lam.ratio.len,ncol=nlambda)

  # for display
  if (lam.ratio.len < 50) display <- FALSE
  if(display){
    cat(paste(c("|",rep(".",50),"|\n","|"),collapse=""))
    display.count <- 1
  }

  # loops
  for (i in (1:lam.ratio.len)){
    #
    lam1 <- lam.ratios[i]
    lam2 <- 1-lam1

    # glmnet fit
    if (bs=="tr"||bs=="rb"){
      fit <- glmnet::glmnet(x=X,y=y,family="poisson",offset=off.const,alpha=alpha,nlambda=nlambda,penalty.factor=c(rep(0,p),rep(lam1,n1-p),rep(lam2,n2)),lambda.min.ratio=glmnet.lam.min)
    } else if (bs=="bs"){
      fit <- glmnet::glmnet(x=X,y=y,family="poisson",offset=off.const,alpha=alpha,nlambda=nlambda,penalty.factor=c(rep(lam1,n1),rep(lam2,n2)),lambda.min.ratio=glmnet.lam.min)
    }

    temp <- try(coef(fit), silent=TRUE)

    if (class(temp)=="try-error"){
      mdl[i,] <- Inf
      mdl.pen[i,] <- Inf
    } else {
      # computing mdl
      nonzeros <- apply(temp,2, function(x,n1,n2){c(sum(x[1:(n1+1)]!=0),sum(x[(n1+2):(n1+n2+1)]!=0))} ,n1=n1,n2=n2)
      if (is.null(assign.emiss)) nonzeros[2,] <- 0
      dev <- (1-fit$dev.ratio)*fit$nulldev
      #mdl.pen[i,] <- log(n*reps)*(fit$df+1)+2*bgamma*lchoose(n2,nonzeros[2,])
      #mdl[i,] <- dev + mdl.pen[i,]
      n3 <- length(fit$lambda)
      lams[i,1:n3] <- fit$lambda
      mdl.pen[i,1:n3] <- log(n*reps)*(fit$df+1)/2+bgamma*lchoose(n2,nonzeros[2,])*emdl
      # since nonzeros[2,] is assigned to 0 if assign.emiss=NULL, the last term will be automatically zero
      mdl[i,1:n3] <- dev/2 + mdl.pen[i,1:n3]*v
      beta.cube[i,1:n3,] <- t(as.matrix(coef(fit)))
    }
    # for display
    if (display){
      if ((i/lam.ratio.len*50)>=display.count) {
        cat(".")
        display.count <- display.count+1
      }
    }
  }
  if (display){
    cat("|\n")
  }

  # output fit
  mdl.ind <- arrayInd(which.min(mdl),dim(mdl))
  if (is.null(assign.emiss)){
    if (mdl.ind[2]==nlambda){
      warnings("mdl: require a larger range\n")
    }
  } else {
    if ((min(mdl.ind)==1)||(mdl.ind[1]==lam.ratio.len)||(mdl.ind[2]==nlambda)){
      warnings("mdl: require a larger range\n")
    }
  }
  best.beta <- as.vector(beta.cube[mdl.ind[1],mdl.ind[2],])
  fitted.logbright <- as.vector(best.beta[1] + X%*%best.beta[-1])
  fitted.logbright.ne <- as.vector(best.beta[1] + B.splines%*%best.beta[2:(n1+1)])
  df <- c(sum(best.beta[2:(n1+1)]!=0),sum(best.beta[(n1+2):(n1+n2+1)]!=0))
  if (is.null(assign.emiss)) df[2] <- 0
  out.ind <- assign.emiss[which(best.beta[(n1+2):(n1+n2+1)]!=0)]
  if (simple){
    beta.cube <- NULL
  }
  return(list(best.beta=best.beta,fitted.logbright=fitted.logbright,fitted.logbright.ne=fitted.logbright.ne,df=df,out.ind=out.ind,mdl.ind=mdl.ind,mdl=mdl,mdl.pen=mdl.pen,beta.cube=beta.cube,reps=reps,off.const=off.const))
}


###########################
#### plotting spectrum ####
###########################
# if there are more utility funciton, we may want to write a Class

plotspec <- function(x, fit, np=TRUE, transform=NULL, ...){
  n <- length(x)
  delta.w <- (x[2]-x[1])
  if (is.null(transform)){
    bright <- exp(fit$fitted.logbright)
    bright.ne <- exp(fit$fitted.logbright.ne)
  } else if (transform=="log"){
    bright <- fit$fitted.logbright
    bright.ne <- fit$fitted.logbright.ne
  } else if (transform=="log10"){
    bright <- fit$fitted.logbright / log(10)
    bright.ne <- fit$fitted.logbright.ne / log(10)
  }
  if (np){
    plot(x,bright,type="n",xlim=c(x[1]-delta.w/2,x[n]+delta.w/2),...)
  }

  # continuum
  lines(x, bright.ne, ...)

  # emission lines
  for (i in (fit$out.ind)){
    if (i==1){
      x2 <- x[1]+delta.w/2
      y2 <- (bright.ne[1]+bright.ne[2])/2
      lines(c(x[1],x2),c(bright[1],bright[1]),...)
      lines(c(x2,x2),c(bright[1],y2),...)
      lines(c(x[1],x2),c(bright.ne[1],y2),...)
    } else if (i==n){
      x1 <- x[n]-delta.w/2
      y1 <- (bright.ne[n-1]+bright.ne[n])/2
      lines(c(x1,x1),c(y1,bright[n]),...)
      lines(c(x1,x[n]),c(bright[n],bright[n]),...)
      lines(c(x1,x[n]),c(y1,bright.ne[n]),...)
    } else {
      x1 <- x[i]-delta.w/2
      x2 <- x[i]+delta.w/2
      y1 <- (bright.ne[i-1]+bright.ne[i])/2
      y2 <- (bright.ne[i]+bright.ne[i+1])/2
      lines(c(x1,x1),c(y1,bright[i]),...)
      lines(c(x1,x2),c(bright[i],bright[i]),...)
      lines(c(x2,x2),c(bright[i],y2),...)
      lines(c(x1,x[i]),c(y1,bright.ne[i]),...)
      lines(c(x[i],x2),c(bright.ne[i],y2),...)
    }
  }
  invisible()
}

