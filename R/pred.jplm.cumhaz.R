pred.jplm.cumhaz <- function(object, at=NULL, CI=FALSE){ 
  digits = max(4, getOption("digits") - 1)
  testpt=ifelse(matrix(object$wt,length(at),length(object$wt),byrow=TRUE) <= matrix(at,length(at),length(object$wt)), 1, 0)
  coef = c(testpt %*% object$lambdas)
  coefTab = cbind(At=at, Value=coef)
  out <- as.data.frame(round(coefTab, digits))
  if (CI==TRUE) { 
  where= length(object$coef.fixed.surv)+length(object$coef.frailty.surv) + (1:length(object$wt)) 
  cov.lambdas = object$covt[where, where]
  se=sqrt(diag((as.matrix(testpt%*%cov.lambdas%*%t(testpt)))))
  lower = coef - 1.96*se
  upper = coef + 1.96*se
  coefTab = cbind(coefTab, Std.Err=se, Lower=lower,Upper=upper, deparse.level = 0)  
  out <- as.data.frame(coefTab)
  ind <- out$Lower < 0
  if(sum(ind)>0)   out$Lower[ind] <- sprintf("<0")
  ind <- out$Upper < 0
  if(sum(ind)>0)   out$Upper[ind] <- sprintf("<0")  
  }
  cat("\n<< Event Process >>\n")
  transf.par=object$K
  if (transf.par==0) cat("Transformation function: H(x) = x\n") else if (transf.par>0)
    cat(sprintf("\nTransformation function: H(x) = log(1+%s)/%s\n",transf.par,transf.par))
  cat("\nBaseline Cumulative Hazard:\n")  
  print(out)  
}
