covest <- function(par_Y, par_T, par_b, K, dl1, dl2, npar, GQl){
  if (K==0){  
  H <- function(x) (x) 
  Hinv <- function(x) {x}   
  H1d <- function(x) 1
  H2d <- function(x) 0
  H3d <- function(x) 0
} else if (K > 0){
  H <- function(x) {log(1+K*x)/K}
  Hinv <- function(x) {(exp(K*x)-1)/K} 
  H1d <- function(x) {1/(1+K*x)}
  H2d <- function(x) {-K/(1+K*x)^2}
  H3d <- function(x) {2*K^2/(1+K*x)^3}
}
   N=npar$N;  nb=length(par_b);
   thetay=par_Y[-npar$nY];  sigma2e= par_Y[npar$nY];    
   gamma=par_T[1:npar$nG];  phi = par_T[npar$nG+(1:npar$nP)];  
   lambda= par_T[(npar$nG+npar$nP)+1:npar$mt];       
   r1= par_b[1]; r2= par_b[2]; r3= par_b[3]; 
   bbv <- Tbb(par_Y,par_T,par_b, K, dl1,dl2, npar, GQl)   
   b1= bbv[1:N, ];              b2=bbv[N+(1:N), ];    
   wwN <- GetWeight(bbv, par_Y, par_T, par_b, K, dl1, dl2, npar, GQl)  
   YAXX1B <- c(dl1$yt- dl1$Xy%*%thetay)  
   L_Vi <- c(dl2$ind_L %*% lambda)     
   exptp <- exp(c(dl2$X2%*%gamma) + phi*b1)	
   expL <- exptp*L_Vi  
 
  H1=H1d(expL);      H2=H2d(expL);       H3=H3d(expL);
  gH2H1 = H2/H1*exptp 
  gdH2H1= H3/H1*exptp*exptp-(H2/H1*exptp)^2 
  gH1 = H1*exptp 	  
  gH2 = H2*exptp*exptp  
  Eb1 <- rowSums(b1*wwN)   
  Eb2 <- rowSums(b2*wwN)   
  Eb1sq <- rowSums(b1^2*wwN) 
  Eb2sq <- rowSums(b2^2*wwN) 
  Eb1b2 <- rowSums(b1*b2*wwN) 
  EH2H1 <- rowSums(gH2H1*wwN)
  EdH2H1 <- rowSums(gdH2H1*wwN)  
  EH1 <- rowSums(gH1*wwN)
  EH2 <- rowSums(gH2*wwN)	
  EH1b1 <- rowSums(gH1*b1*wwN)
  EH2b1<- rowSums(gH2*b1*wwN) 
  EH1b1sq <- rowSums(gH1*b1^2*wwN)
  EH2b1sq <- rowSums(gH2*b1^2*wwN) 
  EH2H1b1 <- rowSums(gH2H1*b1*wwN)
  EdH2H1b1 <- rowSums(gdH2H1*b1*wwN)  
  EH2H1b1sq <- rowSums(gH2H1*b1^2*wwN)
  EdH2H1b1sq <- rowSums(gdH2H1*b1^2*wwN)  

  ED2LY <- matrix(0, npar$nY, npar$nY)    
  ED2LT <- matrix(0, npar$nTl, npar$nTl)       
  ED2Lb <- matrix(0, nb, nb)
  E_DLDL <- matrix(0, (npar$nY+npar$nTl+nb), (npar$nY+npar$nTl+nb))
  EDLEDL <- matrix(0, (npar$nY+npar$nTl+nb), (npar$nY+npar$nTl+nb))    
  for (ii in 1:N){
       IDXY <- which(dl1$id.vec.y==dl2$id.vec[ii])
     D1thetay <- matrix(0, length(thetay), GQl$Nnode)       
     D1sigma2e <- matrix(0, 1, GQl$Nnode)
     D2thetay <- matrix(0, length(thetay), length(thetay)) 
     D2sigma2e <- 0
     D2thetaysigma2e <-matrix(0, length(thetay), 1)           
     for (k in IDXY){
       ### Y
       D1thetay=D1thetay + dl1$Xy[k, ]%*%t((YAXX1B[k]-dl1$ZZ1[k,1]*b1[ii,]-dl1$ZZ1[k,2]*b2[ii,])/sigma2e)       
       D1sigma2e=D1sigma2e + t((YAXX1B[k]-dl1$ZZ1[k,1]*b1[ii,]-dl1$ZZ1[k,2]*b2[ii,])^2/(2*sigma2e^2)-1/(2*sigma2e))
       D2thetay <- D2thetay -dl1$Xy[k, ]%*%t(dl1$Xy[k, ])/sigma2e    
       D2thetaysigma2e=D2thetaysigma2e -(YAXX1B[k]-dl1$ZZ1[k,1]*Eb1[ii]-dl1$ZZ1[k,2]*Eb2[ii])/sigma2e^2*dl1$Xy[k, ]        
       temp=(YAXX1B[k]-dl1$ZZ1[k,1]*b1[ii,]-dl1$ZZ1[k,2]*b2[ii,])^2
       D2sigma2e=D2sigma2e + 1/2/sigma2e^2-sum(temp*wwN[ii,])/sigma2e^3      
     }     
     RISK <- which(dl2$wt <= dl2$EventTime[ii]) 
     indX <- ifelse(dl2$wt <= dl2$EventTime[ii], 1, 0)      
     coef2 <- dl2$DELTA[ii]*gH2H1[ii, ]-gH1[ii, ]       
     coef4 <- dl2$DELTA[ii]*(EH2H1[ii]+EdH2H1[ii]*L_Vi[ii])-(EH1[ii]+EH2[ii]*L_Vi[ii])     
     coef5b1=dl2$DELTA[ii]*(EH2H1b1[ii]+EdH2H1b1[ii]*L_Vi[ii])-(EH1b1[ii]+EH2b1[ii]*L_Vi[ii])                                    
     D1G <-  dl2$DELTA[ii]*dl2$X2[ii, ] + dl2$X2[ii, ]%*%t(coef2*L_Vi[ii])
     D1phi1 <- t(dl2$DELTA[ii]*b1[ii,] + coef2*b1[ii,]*L_Vi[ii])         
     D1lambda <- dl2$DELTA[ii]*ifelse(dl2$wt==dl2$EventTime[ii],1,0)/lambda + indX%*%t(coef2)
     D2G <- coef4*L_Vi[ii]*dl2$X2[ii, ]%*%t(dl2$X2[ii, ])              
     D2Gphi1 <- matrix(coef5b1*L_Vi[ii]*dl2$X2[ii, ])         
     D2Glambda <- matrix(0, npar$nG, npar$mt)
     D2Glambda[ ,RISK] <- coef4*dl2$X2[ii, ]
     D2phi1 <- -(EH1b1sq[ii]+EH2b1sq[ii]*L_Vi[ii])*L_Vi[ii]
         if (dl2$DELTA[ii]==1) {D2phi1 =D2phi1+(EH2H1b1sq[ii]+EdH2H1b1sq[ii]*L_Vi[ii])*L_Vi[ii]} 
     D2phi1lambda <- matrix(0, 1, npar$mt)
     D2phi1lambda[ ,RISK] <- coef5b1  
     D2lambda <- indX%*%t(indX)*(dl2$DELTA[ii]*EdH2H1[ii]-EH2[ii]) 
     if (dl2$DELTA[ii]==1) {D2lambda = D2lambda-diag(npar$mt)*ifelse(dl2$wt==dl2$EventTime[ii],1,0)/lambda^2}   
     # random effect      
     cb1 <- t((b1[ii,]/r1)^2)
     cb2 <- t((b2[ii,]/r2)^2)      
     cb12 <- t(r3*b1[ii,]*b2[ii,]/(r1*r2)^0.5) 
     Ecb1 <- Eb1sq[ii]/r1^2
     Ecb2 <- Eb2sq[ii]/r2^2
     Ecb1b2 <- r3*Eb1b2[ii]/(r1*r2)^0.5    
     Dr1 <- -1/2/r1 -1/2/(1-r3^2)*(-cb1 + cb12/r1)     
     Dr2 <- -1/2/r2 -1/2/(1-r3^2)*(-cb2 + cb12/r2)     
     Dr3 <- (r3-r3/(1-r3^2)*(r1*cb1-2*cb12+r2*cb2)+cb12/r3)/(1-r3^2)
     D2r1<- 1/2/r1^2 - (2*Ecb1/r1-1.5*Ecb1b2/r1^2)/(2*(1-r3^2))     
     D2r2<- 1/2/r2^2 - (2*Ecb2/r2-1.5*Ecb1b2/r2^2)/(2*(1-r3^2)) 
     D2r3<- ((1+r3^2) - (1+3*r3^2)/(1-r3^2)*(Ecb1*r1-2*Ecb1b2+Ecb2*r2) + 4*Ecb1b2)/(1-r3^2)^2
     D2r1r2<- Ecb1b2/(r1*r2*4*(1-r3^2))
     D2r1r3<- r3/(1-r3^2)^2*(Ecb1-Ecb1b2/r1)-Ecb1b2/r3/(2*(1-r3^2)*r1)
     D2r2r3<- r3/(1-r3^2)^2*(Ecb2-Ecb1b2/r2)-Ecb1b2/r3/(2*(1-r3^2)*r2)               
     
     ED2LY <- ED2LY + rbind(cbind(D2thetay, D2thetaysigma2e), cbind(t(D2thetaysigma2e), D2sigma2e))
     ED2LT <- ED2LT + rbind(cbind(D2G, D2Gphi1, D2Glambda), 
                       	    cbind(t(D2Gphi1), D2phi1, D2phi1lambda),
                       	    cbind(t(D2Glambda), t(D2phi1lambda), D2lambda))
     ED2Lb <- ED2Lb + cbind(c(D2r1,D2r1r2,D2r1r3), c(D2r1r2,D2r2,D2r2r3), c(D2r1r3,D2r2r3,D2r3))	 	                 
     DL <- rbind(D1thetay,D1sigma2e, D1G,D1phi1,D1lambda, Dr1,Dr2,Dr3)
     EDL <-  t(DL)*wwN[ii,]
     E_DLDL <- E_DLDL + DL%*%(EDL)
     EDL <- colSums(EDL) 
     EDLEDL <- EDLEDL + EDL%*%t(EDL) 
   }
   ED2L <- rbind(cbind(ED2LY, matrix(0,npar$nY,npar$nTl), matrix(0,npar$nY,nb)),                
                 cbind(matrix(0,npar$nTl,npar$nY), ED2LT, matrix(0,npar$nTl,nb)),
                 cbind(matrix(0,nb,npar$nY), matrix(0,nb,npar$nTl), ED2Lb))
   HS <- -ED2L - (E_DLDL-EDLEDL)  
   COV <- solve(HS); COV;  
   errCOV=0;   
   if (length(which(diag(COV)<0))>0){    
   out.eigen<-eigen(HS)
   values<-out.eigen$values	
   temp1<-which(values < 0)
   if (length(temp1) > 0){
   		values[-temp1]<-1/values[-temp1]
   		values[temp1]<-0
            errCOV=1;
   	} else {values <- 1/values}
   COV <- (out.eigen$vectors) %*% diag(values) %*% t(out.eigen$vectors);
   }   
   outCOV <- list(COV=COV, errCOV=errCOV);  outCOV;   
}
