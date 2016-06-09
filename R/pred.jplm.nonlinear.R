pred.jplm.nonlinear <- function(object, nlm.par, at=NULL, CI=FALSE){
  digits = max(4, getOption("digits") - 4)	
  hatbasis = basis(at, object$degree,object$n.knots,bdr=c(min(nlm.par), max(nlm.par)))	
  coef = c(hatbasis %*% object$coef.nlm.y)
  cov.zetas = object$covy[length(object$coef.lm.y)+(1:length(object$coef.nlm.y)), length(object$coef.lm.y)+(1:length(object$coef.nlm.y))]
  se=sqrt(diag((as.matrix(hatbasis%*%cov.zetas%*%t(hatbasis)))))
  coefTab = cbind(At=at, Value=coef, Std.Err=se, deparse.level = 0)  
  cat("\n<< Longitudinal Process >>\n")
  cat(sprintf("B-spline estimation was used with degree=%s and %s inner knots.\n", object$degree, object$n.knots))
  cat("\nNonlinear Coefficient:\n")
if (CI){
  lower = coef - 1.96*se
  upper = coef + 1.96*se
  coefTab = cbind(coefTab, Lower=lower, Upper=upper,	deparse.level = 0)
  out <- as.data.frame(coefTab)
} else{
  coefTab = cbind(coefTab, z_value=coef/se,p_value=2*pnorm(abs(coef/se),lower.tail=FALSE), deparse.level = 0) 	
  out <- as.data.frame(round(coefTab, digits)) 
  ind <- out$p_value == 0
  out$p_value <- sprintf(paste("%.", digits, "f", sep = ""), out$p_value)
  out$p_value [ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")	
}  
  print(out)  
}