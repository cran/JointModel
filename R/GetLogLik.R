GetLogLik_phi_b1 <- function(par_Y, par_T, par_b, K, dl1, dl2, npar, GQl){
  if (K==0){  
  H <- function(x) (x) 
  H1d <- function(x) 1
} else if (K > 0){
  H <- function(x) {log(1+K*x)/K}
  H1d <- function(x) {1/(1+K*x)}
}

   N=npar$N
   thetay=par_Y[-npar$nY];  sigma2e= par_Y[npar$nY];    
   gamma=par_T[1:npar$nG];  phi = par_T[npar$nG+(1:npar$nP)];  
   lambda= par_T[(npar$nG+npar$nP)+1:npar$mt]; 
   YAXX1B <- c(dl1$yt- dl1$Xy%*%thetay)   
   Ycmp <- dl1$indIDY%*%(-log(sigma2e)/2 -(YAXX1B)^2/(2*sigma2e)) 

  Sb =  par_b;   ### par_b[3]=corr(b1, b2) == r3
  if (length(par_b) > 1){
	Sb = diag(par_b[1:2])
	Sb[1,2]= par_b[3]*sqrt(par_b[1]*par_b[2]); 
	Sb[2,1]=Sb[1,2]
  }
  eV <- eigen(Sb)
  Sb.sq<- eV$vectors %*% diag(sqrt(eV$values)) %*% t(eV$vectors)
  Sb.inv<-eV$vectors %*% diag(1/eV$values) %*% t(eV$vectors)
  bcmp2 <- rep(0,N)
  for (i in 1:N){
     idx = which(dl1$indIDY[i, ]==1)
     if (length(idx)>1){
        Di.inv = Sb.inv + t(dl1$ZZ1[idx, ])%*%dl1$ZZ1[idx, ]/sigma2e
        mui=t(dl1$ZZ1[idx, ])%*%YAXX1B[idx]/sigma2e + c(dl2$DELTA[i]*phi, 0)      
     } else{
     	Di.inv = Sb.inv + dl1$ZZ1[idx, ]%*%t(dl1$ZZ1[idx, ])/sigma2e
        mui=dl1$ZZ1[idx, ]*YAXX1B[idx]/sigma2e + c(dl2$DELTA[i]*phi, 0) 
     }	
     ev2 = eigen(Di.inv)
     Di =ev2$vectors%*%diag(1/ev2$values)%*%t(ev2$vectors)
     bcmp2[i] <- t(mui) %*% Di %*% mui
  } 
  bcmp <- -log(det(Sb))/2 + 0.5*bcmp2

  bbv <- Tbb(par_Y, par_T, par_b, K, dl1, dl2, npar, GQl)   
  b1= bbv[1:N, ];              b2=bbv[N+(1:N), ];    
  bby = b1[dl1$id.vec.y, ] + b2[dl1$id.vec.y, ]*dl1$ZZ1[,2]
  exptp = exp(c(dl2$X2%*%gamma) + phi*b1)  
  L_Vi = c(dl2$ind_L %*% lambda)
  expL <- exptp*L_Vi      ## N x Nnode matrix
  hi_b=exp(-H(expL))  
  hi_b <- hi_b*(H1d(expL)^(dl2$DELTA)) 
  Ehi_b <- rowSums(hi_b*matrix(GQl$ww.prod, N, GQl$Nnode, byrow=TRUE))  
  Tcmp2 = dl2$DELTA*c(dl2$X2%*%gamma) + log(Ehi_b)  
  loglik <- sum(log(lambda)) + sum(Ycmp + bcmp + Tcmp2); loglik;   
}
