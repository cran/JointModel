
quad <- function(db, NGHQ=4){
  GHQ <- gauss.quad(NGHQ, "hermite")	
  node <- GHQ$nodes
  weight <- GHQ$weights 
 
# dim(b): db==2;
  if (db==2){
     uu = rbind(rep(node, rep(NGHQ,NGHQ)),rep(node,NGHQ))
     ww = rbind(rep(weight, rep(NGHQ,NGHQ)),rep(weight,NGHQ))
     ww.prod = apply(ww, 2, prod)
     Nnode <- dim(ww)[2]
  }
  out <- list(uu=uu, ww.prod=ww.prod, Nnode=Nnode); 
  out;
}

### SK: revise to be used w/ db=1
Tbb <- function(par_Y, par_T, par_b, K, dl1, dl2, npar, GQl){
  thetay=par_Y[-npar$nY];  sigma2e= par_Y[npar$nY];  
  phi = par_T[npar$nG+(1:npar$nP)];    
  Sb =  par_b;   ### par_b[3]=corr(b1, b2) == r3
  if (npar$db > 1){
	Sb = diag(par_b[1:npar$db])
	Sb[1,2]= par_b[3]*sqrt(par_b[1]*par_b[2]); 
	Sb[2,1]=Sb[1,2]
  }
  eV <- eigen(Sb)
  Sb.inv<- eV$vectors %*% diag(1/eV$values) %*% t(eV$vectors)
  YAXX1B <- c(dl1$yt- dl1$Xy%*%thetay) 
    
  bbv <- matrix(0, npar$N*npar$db, GQl$Nnode)    
  for (i in 1:npar$N){
     idx = which(dl1$indIDY[i, ]==1)
     if (length(idx)>1){
        Vb.inv = Sb.inv + t(dl1$ZZ1[idx, ])%*%dl1$ZZ1[idx, ]/sigma2e
        mu = t(dl1$ZZ1[idx, ])%*%YAXX1B[idx]/sigma2e + c(dl2$DELTA[i]*phi, 0)
     } else{
     	Vb.inv = Sb.inv + dl1$ZZ1[idx, ]%*%t(dl1$ZZ1[idx, ])/sigma2e
     	mu = dl1$ZZ1[idx, ]*YAXX1B[idx]/sigma2e + c(dl2$DELTA[i]*phi, 0)
     }	
     ev2 = eigen(Vb.inv)
     Vb =ev2$vectors%*%diag(1/ev2$values)%*%t(ev2$vectors)
     Vb.sq =ev2$vectors%*%diag(1/sqrt(ev2$values))%*%t(ev2$vectors)     
     
     bb = sqrt(2)*Vb.sq%*%GQl$uu + matrix(Vb%*%mu, npar$db, GQl$Nnode)
     bbv[c(i, npar$N+i), ] = bb 
  }; bbv;  
}

GetWeight <- function(bbv, par_Y, par_T, par_b, K, dl1, dl2, npar, GQ){
  if (K==0){  
  H <- function(x) {x} 
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
  N = length(dl2$DELTA)	
  gamma =par_T[1:npar$nG];  phi = par_T[npar$nG+(1:npar$nP)];  
  lambda = par_T[(npar$nG+npar$nP)+1:npar$mt];    
  b1= bbv[1:N, ];     b2=bbv[N+(1:N), ];
  exptp <- exp(c(dl2$X2%*%gamma) + phi*b1) 
  L_Vi <- c(dl2$ind_L %*% lambda)  
  out=exp(-H(exptp*L_Vi))  
  out <- out*(H1d(exptp*L_Vi)^(dl2$DELTA))      
  weights<- out*matrix(GQ$ww.prod, N, GQ$Nnode, byrow=TRUE)
  wwN <- weights/ rowSums(weights);  wwN; 
} 

M_step <- function(bbv, wwN, par_Y, par_T, par_b, K, dl1, dl2, npar){ 
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

   N=npar$N
   thetay=par_Y[-npar$nY];  sigma2e= par_Y[npar$nY];       
   gamma=par_T[1:npar$nG];   phi = par_T[npar$nG+(1:npar$nP)];  
   lambda= par_T[(npar$nG+npar$nP)+1:npar$mt];  
   YAXX1B <- c(dl1$yt- dl1$Xy%*%thetay) 
   L_Vi <- c(dl2$ind_L %*% lambda)  
   b1= bbv[1:N, ];     b2=bbv[N+(1:N), ];
   exptp <- exp(c(dl2$X2%*%gamma) + phi*b1)
   q <- -dl2$DELTA*H2d(exptp*L_Vi)/H1d(exptp*L_Vi) + H1d(exptp*L_Vi)      

   ####### Expectations
   Eb1 <- rowSums(b1*wwN)   
   Eb2 <- rowSums(b2*wwN)   
   Eb1sq <- rowSums(b1^2*wwN) 
   Eb2sq <- rowSums(b2^2*wwN) 
   Eb1b2 <- rowSums(b1*b2*wwN) 
   Eqe <-  rowSums(q*exptp*wwN)	
   Eqeb1 <-rowSums(q*exptp*b1*wwN)
   Eqeb1sq <-rowSums(q*exptp*b1^2*wwN)
   
   ####### Update param set 1
   newthetay=solve(t(dl1$Xy)%*%dl1$Xy)%*%t(dl1$Xy)%*%(dl1$yt-Eb1[dl1$id.vec.y] - Eb2[dl1$id.vec.y]*dl1$ZZ1[,2])
   temp.e = YAXX1B - b1[dl1$id.vec.y, ] - b2[dl1$id.vec.y, ]*dl1$ZZ1[,2]
   newsigma2_e= sum(temp.e^2 *wwN[dl1$id.vec.y, ]) /length(YAXX1B)
   newr1=sum(Eb1sq)/N 
   newr2=sum(Eb2sq)/N 
   newr3=sum(Eb1b2)/N /sqrt(newr1*newr2)

   ####### Update param set 2
   sEqe <- c(dl2$ind_VjVi%*%Eqe) 
   sEqeb1 <- c(dl2$ind_VjVi%*%Eqeb1)/sEqe
   sXEqe <- (dl2$ind_VjVi%*%(dl2$X2*Eqe))/sEqe 
   g1 <- colSums(dl2$DELTA*(dl2$X2 - sXEqe))
   g2 <- sum(dl2$DELTA*(Eb1 - sEqeb1))
   XEqe <- (dl2$X2*Eqe)/sEqe                  
   tempG<- array(0, c(npar$nG,npar$nG,N))  
   dg1dG<- matrix(0, npar$nG, npar$nG)
   for (ii in 1:N){
     tempG[,,ii]<-dl2$X2[ii, ] %*% t(XEqe[ii, ])  
   } 
   idxD <- which(dl2$DELTA==1)
   for (kk in idxD){
     RISK <- which(dl2$EventTime>=dl2$EventTime[kk])
     comp1G <-0
     for (jj in RISK){
	      comp1G <- comp1G + tempG[,,jj]
     }
     dg1dG <- dg1dG - (comp1G - sXEqe[kk, ]%*%t(sXEqe[kk, ]))
   }    
   dg1dP <- colSums(-dl2$DELTA*((dl2$ind_VjVi%*%(dl2$X2*Eqeb1))/sEqe - sXEqe*sEqeb1))    
   dg2dG <- t(dg1dP)      
   dg2dP <- sum(-dl2$DELTA*((dl2$ind_VjVi%*%Eqeb1sq)/sEqe - sEqeb1^2))    
   HH <- rbind(  cbind(dg1dG,dg1dP),
                 cbind(dg2dG,dg2dP)) 
   param <- c(gamma, phi)
   newparam <- param - solve(HH)%*%c(g1, g2) 
   
   ####### Update param set 3 - lambda
   newlambda = dl2$d_wj / c(t(dl2$ind_L) %*% Eqe)

   newpar_Y <- c(newthetay, newsigma2_e)
   newpar_T <- c(newparam, newlambda)
   newpar_b <- c(newr1, newr2, newr3)    
   c(newpar_Y, newpar_T, newpar_b)   
}   # end of "M_step" # 

