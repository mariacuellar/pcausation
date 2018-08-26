
#' Semiparametric estimation of the probability of causation
#'
#' Takes in values of binary outcome y, binary exposure a, and covariates x1...xn, and outputs prob of causation and other parameters.
#' @param x A data frame of covariates x1...xn.
#' @param y A binary vector, of 0s and 1s, representing the outcome (e.g. of an intervention).
#' @param a A binary vector, of 0s and 1s, representing the exposure (e.g. to a treatment).
#' @param xtest A data frame of covariates x1...xn at which to evaluate the probability of causation. Could be same as x.
#' @param nsplits How many splits to use for cross-validation. If none are desired, use 1.
#' @param start.list Which values to use in starting optimization. Default is vector of 0s.
#' @param tracetf Do you want to show the trace, i.e., progress, of the optimization?
#' @param printres Do you want to print the table of results?
#' @param showprogress Do you want to see a progress bar while the estimation is performed?
#' @param sl.lib Library of algorithms to be given to superlearner. Default is a mixture.
#' @return A table of estimates of the probability of causation with confidence intervals, and other parameters.
#' @export


pcausation = function(y, a, x, xtest, nsplits=2, start.list=c(rep(0,ncol(x))), tracetf=FALSE, printres=TRUE, showprogress=TRUE,
                      sl.lib = c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.rpart")){

  # packages
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("rpart")

  # setting things up
  n <- nrow(x)
  avals <- names(table(a))
  n.avals <- length(avals)
  s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

  # progress bar
  if(showprogress==TRUE){ pb <- txtProgressBar(min=0, max=2*nsplits*n.avals, style=3) }

  muhat <- as.data.frame(matrix(NA,nrow=n,ncol=n.avals))
  colnames(muhat) <- paste("a",avals,sep="")
  pihat <- muhat
  muhat.xtest <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=n.avals))
  pihat.xtest <- muhat.xtest

  # estimating nuisance parameters: propensity score and outcome regressions
  pbcount <- 0
  for (i in 1:n.avals){
    if (i==1){ Sys.sleep(0.1); if(showprogress==TRUE){setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }}

    # starting cross validation
    for (vfold in 1:nsplits){

      # split into train and test sets
      train <- s!=vfold; test <- s==vfold
      if (nsplits==1){ train <- test }

      # estimate propensity score (pi)
      if (i != n.avals){

        pifit.xtest <- SuperLearner(as.numeric(a==avals[i])[train],as.data.frame(x[train,]),
                                    newX=xtest, SL.library=sl.lib, family=binomial)

        pihat.xtest <- as.numeric(pifit.xtest$SL.predict) # evaluated at pre-selected set of x's

        Sys.sleep(0.1)
        if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
      }

    }

    # estimate outcome regression function (mu), one for each value of a
    mufit.xtest <- SuperLearner(y[a==avals[i] & train],
                                as.data.frame(x[a==avals[i] & train,]),
                                newX=xtest, SL.library=sl.lib)

    muhat.xtest[,i] <- mufit.xtest$SL.predict # evaluated at pre-selected set of x's

    Sys.sleep(0.1)
    if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }

  }
  if (i == n.avals){ pihat[,i] <- 1 - apply(pihat,1,sum, na.rm=T) }

  # estimates of E{Y(0)} and E{Y(1)}
  names(muhat.xtest) <- c("a0", "a1")
  amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
  alevel <- matrix(rep(as.numeric(avals), rep(n,n.avals)),nrow=n,byrow=F)
  ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)

  # estimates of probability of causation - plugin estimation
  est.pi <- abs(1 - muhat.xtest[1]/muhat.xtest[2]); names(est.pi)<-"est.pi" # se and CI's are not of interest

  # estimates of probability of causation - influence function estimation
  include.x = names(as.data.frame(x))
  include.b = gsub("x", "b", include.x)
  include.xb = diag(sapply( include.b, FUN=function(x, y) paste(x, include.x, sep="*") ))
  formula.x = paste("y ~ ", paste(include.x, collapse="+"), sep="")
  formula.pi = paste("a ~ ", paste(include.x, collapse="+"), sep="")
  formula.nls.x = paste("y ~ expit(", paste(include.xb, collapse="+"), ")", sep = "")
  start.values = setNames(as.list(start.list), c(include.b))
  ystar = (1/muhat.xtest[2])*((muhat.xtest[1]/muhat.xtest[2])*(1/pihat.xtest)*a*(y-muhat.xtest[2]) -
                                (1/(1-pihat.xtest))*(1-a)*(y-muhat.xtest[1])) + (1 - muhat.xtest[1]/muhat.xtest[2])
  ifvals <- ystar; names(ifvals) = "ystar"

  # fit nls model - gives NA if the NLS doesn't converge
  mod = tryCatch({
    nls( as.formula(formula.nls.x), start=start.values, data=as.data.frame(cbind(x,y,a)), nls.control(maxiter = 500), algorithm = "port",
         lower=c(rep(-5, length(include.x))), upper=c(rep(5, length(include.x))), trace = tracetf )
  }, warning = function(warning_condition) {
    mod <- NA
  }, error = function(error_condition) {
    mod <- NA
  })

  # estimate covariance matrix for betahat
  preds <- try( predict(mod) , silent = TRUE )
  xmat <- try( as.matrix(x) , silent = TRUE ); try( xtestmat <- as.matrix(xtest) , silent = TRUE )
  wts <- try( mod$weights , silent = TRUE ); try( if (is.null(wts)){ wts <- 1 } , silent = TRUE )
  bread <- try( solve( (t(xmat *(preds*(1-preds)*wts)) %*% xmat)/n ) , silent = TRUE )
  meat <- try( (t(xmat *( ((y-preds) * wts)^2)) %*% xmat)/n , silent = TRUE )
  vcov <- try( bread %*% meat %*% bread / n , silent = TRUE )
  coefs <- try( coef(mod) , silent = TRUE )
  res.betas <- tryCatch({ data.frame(Estimate=coefs, Robust.SE=sqrt(diag(vcov)),
                                     z.val=coefs/sqrt(diag(vcov)),	p.val= round(2*(1-pnorm(abs(coefs/sqrt(diag(vcov))))),3) ,
                                     ci.ll= coefs-1.96*sqrt(diag(vcov)) , ci.ul=coefs+1.96*sqrt(diag(vcov)) )
  }, warning = function(warning_condition){ res.betas })

  # get predicted value / CI for probability at specific x0
  x0 <- xtestmat
  est.if <- try( expit(x0 %*% coefs) , silent = TRUE )
  se.if <- try( as.data.frame(apply(  x0, 1, function(M){(t(M) %*% vcov %*% M)}  )) , silent = TRUE ); names(se.if)="se.if"
  ci.ll.if = try( as.data.frame(apply(  x0, 1, function(M){expit(M %*% coefs - 1.96*(t(M) %*% vcov %*% M))}  )) , silent = TRUE ); names(ci.ll.if)="ci.ll.if"
  ci.ul.if = try( as.data.frame(apply(  x0, 1, function(M){expit(M %*% coefs + 1.96*(t(M) %*% vcov %*% M))}  )) , silent = TRUE ); names(ci.ul.if)="ci.ul.if"

  # get SEs/CIs/pvals for betas
  coefs <- try( coef(mod), silent=TRUE)
  res.coefs <- try( data.frame(Estimate=coefs, Robust.SE=sqrt(diag(vcov)),
                               z.val=coefs/sqrt(diag(vcov)),	p.val= round(2*(1-pnorm(abs(coefs/sqrt(diag(vcov))))),3) ,
                               ci.ll= coefs-1.96*sqrt(diag(vcov)) , ci.ul=coefs+1.96*sqrt(diag(vcov)) ), silent=TRUE)

  # table of values for influence function estimator and plugin estimator
  res.if <- tryCatch({
    data.frame( est.pc = est.if, se = se.if, ci.ll = ci.ll.if, ci.ul = ci.ul.if )
  }, warning = function(warning_condition) {
    empty = rep(NA,nrow(dat.eval))
    data.frame( est.pc = empty, se = empty, ci.ll = empty, ci.ul = empty )
  }, error = function(error_condition) {
    empty = rep(NA,nrow(dat.eval))
    data.frame( est.pc = empty, se = empty, ci.ll = empty, ci.ul = empty )
  })

  res.pi <- data.frame( est.pc = est.pi )

  Sys.sleep(0.1)
  if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount) ; close(pb) }

  nuis <- as.data.frame(cbind(pihat.xtest,muhat.xtest))
  colnames(nuis) <- c("pihat", "mu0hat", "mu1hat")

  if(printres==TRUE){ print(res.if) }
  return(invisible(list(res.if=res.if, res.pi=res.pi, nuis=nuis, ifvals=ifvals, res.coefs=res.coefs)))

}
