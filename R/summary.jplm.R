summary.jplm <- function(object, digits=max(4,getOption("digits")-4)){
	if (!inherits(object, "jplm"))
	   stop("Use only with 'jplm' objects \n")
  
  betas <- object$coef.lm.y
  zetas <- object$coef.nlm.y
  gammas <- object$coef.fixed.surv
  phis <- object$coef.frailty.surv
  sigma2_e <-  object$var.resid
  ran.vcomp <- object$raneff.vcomp

  if( is.null(object$covy)){
  tab.y <- data.frame("Value"=betas, row.names=object$vnf.y)
  tab.t1 <- as.data.frame(cbind("Value"=gammas, "exp(Value)"=exp(gammas)), row.names=object$vnf.t)  
  tab.t2 <- as.data.frame(cbind("Value"=phis, "exp(Value)"=exp(phis)), row.names='phi') 
  } else{
  sd.betas = sqrt(diag(as.matrix(object$covy[1:length(betas),1:length(betas)])))
  sd.gammas= sqrt(diag(as.matrix(object$covt[1:length(gammas),1:length(gammas)])))
  where = length(gammas) + (1:length(phis))
  sd.phis = sqrt(diag(as.matrix(object$covt[where,where])))	
  coefTab.y <- cbind("Value"=betas, "Std.Err"=sd.betas,"z-value"=betas/sd.betas, 
          "p-value"=2*pnorm(abs(betas/sd.betas),lower.tail=FALSE))
  rownames(coefTab.y) <- object$vnf.y
  coefTab.t1 <- cbind("Value"=gammas, "exp(Value)"=exp(gammas), "Std.Err"=sd.gammas, 
         "p-value"=2*pnorm(abs(gammas/sd.gammas), lower.tail=FALSE))
  rownames(coefTab.t1) <- object$vnf.t        
  coefTab.t2 <- cbind("Value"=phis, "exp(Value)"=exp(phis), "Std.Err"=sd.phis, 
        "z-value"=phis/sd.phis, "p-value"=2*pnorm(abs(phis/sd.phis), lower.tail=FALSE)) 
  rownames(coefTab.t2) <- 'phi'        
  tab.y <- as.data.frame(round(coefTab.y, digits)) 
  ind1 <- tab.y$"p-value" == 0
  tab.y$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), tab.y$"p-value")
  tab.y$"p-value"[ind1] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
  tab.t1 <- as.data.frame(round(coefTab.t1, digits)) 
  ind2 <- tab.t1$"p-value" == 0
  tab.t1$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), tab.t1$"p-value")
  tab.t1$"p-value"[ind2] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "") 
  tab.t2 <- as.data.frame(round(coefTab.t2, digits)) 
  ind2 <- tab.t2$"p-value" == 0
  tab.t2$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), tab.t2$"p-value")
  tab.t2$"p-value"[ind2] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
  }
	
  D <- object$covb
  db <- length(object$vnr.y)
  coef = round(c(sigma2_e, ran.vcomp), digits) 
  attr1 = rep("Variance", 1+db)
  namrow <- c("Residual", object$vnr.y)
  if (db>1){
      attr1 = c(attr1, rep("Corr", length(coef)-(1+db)))	
      namcorr <- object$vnr.y[1]	
      for (i in 2:db){
	    namcorr  <- paste(namcorr ,", ",object$vnr.y[i],sep="")
      }
      namrow <- c(namrow, namcorr)
  }
  tab.vcomp = data.frame("Attr"=attr1, "Value"=coef,  row.names= namrow)
  model.info = data.frame(loglik=object$loglik, AIC=object$AIC, BIC=object$BIC, row.names="Value")
 
  cat("\n<< Longitudinal Process >>\n"); 
  cat("Formula: "); print(object$formula.y); 
  cat("\nLinear Effects Coefficients:\n")
  print(tab.y)
  cat("\nRandom Effects Variance Components:\n")
  print(tab.vcomp)
  cat("\n<< Event Process >>\n"); 
  transf.par=object$K
  if (transf.par==1) cat("Transformation function: H(x) = log(1 + x)\n") else if (transf.par==0)
     cat("Transformation function: H(x) = x\n") else if (transf.par>0)
     cat(sprintf("Transformation function: H(x) = log(1+%s*x)/%s\n",transf.par,transf.par))
  cat("\nProportional Hazards Coefficients (fixed effects):\n")
  print(tab.t1)
  cat("\nProportional Hazards Coefficients (random effects):\n")
  print(tab.t2)
  cat('\n')
  print(model.info)
}