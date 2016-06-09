jplm <- 
function(formula.lm.y, nlm.par=NULL, data.y, formula.surv.fixed, formula.frailty, id.vec=NULL, transf.par=0, data.surv, 
degree=3, n.knots=3, CovEst=TRUE, maxiter = 200, epsilon=5e-04,...){
	
#######
####### survival
#######
if (!is.vector(id.vec)){
	stop("\nError: argument 'id.vec' must be a vector.")
}
if (transf.par<0){
	stop("\nError: argument 'transf.par' must be a non-negative value.")
}

fit.coxph <- coxph(formula.surv.fixed, data = data.surv)
obj2 <- fit.coxph
tt <- delete.response(terms(obj2))
mf2 <- model.frame(obj2, data=data.surv)
temp = model.extract(frame=mf2, "response");
EventTime = temp[,1]
DELTA = temp[,2]
N = length(EventTime)
X2 = as.matrix(model.matrix(obj2, data=data.surv))
Z2 <- as.matrix(model.matrix(formula.frailty, data=data.surv))
varname.fixed.t <- colnames(X2)	
dlt <- list(id.vec=id.vec, EventTime=EventTime,DELTA=DELTA,X2=X2,Z2=Z2)

#######
####### longitudinal
#######
obj1<-lmer(formula=formula.lm.y, data=data.y, REML=FALSE)
yt <- getME(obj1, "y")
XX1 = getME(obj1, "X")     ### Or, XX1=model.matrix(fit.lmer) 
cname=colnames(XX1)
tt <- delete.response(terms(obj1,fixed.only=TRUE))
if(attr(tt, "intercept")==1){
	XX1=as.matrix(XX1[,-1])
	colnames(XX1) <- cname[-1]
}
varname.fixed.y=attr(tt, "term.label")

tt <- delete.response(terms(obj1,random.only=TRUE))
rfd <- model.frame(tt, data.y)
id.vec.y <- rfd[,ncol(rfd)]
if(attr(tt, "intercept")==1){
	ZZ1=matrix(1,nrow(rfd),1); namcol='(Intercept)'
} else {ZZ1=NULL; namcol=NULL}
if (ncol(rfd)>1){
	ZZ1=cbind(ZZ1, rfd[,-ncol(rfd)])
	namcol=c(namcol, colnames(rfd)[1:(ncol(rfd)-1)])
}
colnames(ZZ1) <- namcol
db=ncol(ZZ1)
if (db!=2){
	stop("\nError: current version only allows two-dimensional random effects.")
} 
varname.random.y <- colnames(ZZ1)

### Construct a sieve space of nonlinear trajectory
basis.matrix=basis(nlm.par,degree,n.knots,bdr=c(min(nlm.par), max(nlm.par)))
n.basis = ncol(basis.matrix)
Xy = cbind(XX1, basis.matrix)
dly <- list(id.vec.y=id.vec.y,yt=yt,Xy=Xy,ZZ1=ZZ1)

### Set other algorithm variables
ni <- rep(0, N)
for(i in 1:N){
		ni[i] <- length(which(id.vec.y == id.vec[i]))
} 
dly$ni = ni;   
dlt$wt <- unique(sort(EventTime[which(DELTA==1)]));   mt <- length(dlt$wt)     	 	
d_wj = rep(0, length(dlt$wt))
for (j in 1:mt){
    idx.wt=which(dlt$wt[j]==EventTime*DELTA)
    d_wj[j] = length(idx.wt)
}
dlt$d_wj =d_wj;  ## don't need in jplm

indIDY <- matrix(0, N, length(yt))   ## check: where should I put this func?
for (i in 1:N){
	indIDY[i, which(id.vec.y==id.vec[i])] <- 1
}
dly$indIDY <- indIDY; ## don't need in jplm
dlt$ind_L <- ifelse(matrix(dlt$wt,N,mt,byrow=TRUE) <= matrix(EventTime,N,mt), 1, 0) 
dlt$ind_VjVi <- ifelse(matrix(EventTime,N,N,byrow=TRUE) >= matrix(EventTime,N,N), 1, 0)

# # of params
nB = ncol(XX1);   nG = ncol(X2);   nP=1
nY= nB + n.basis + 1;
nb=db*(db+1)/2;  
nTl = nG+nP+mt; # # of params for covariance matrix
npar <- list(nB=nB, nY=nY, nG=nG, nP=nP, nTl=nTl, db=db, N=N, mt=mt)

GQ <- quad(db=db, NGHQ=4)

###	Initials
thetay0 <- solve(t(Xy)%*%Xy)%*%t(Xy)%*%yt
vcomp=as.data.frame(VarCorr(obj1))$vcov
sigma2_e0=1;
gamma0=obj2$coef;  phi0= rep(0, nP);  lambda0=rep(1/mt, mt);
par_Y <- c(thetay0, sigma2_e0)
par_T <- c(gamma0, phi0, lambda0)
par_b <- rep(0, nb);   par_b[1:db] <- vcomp[1:db]

n.conv <- NULL
for (iter in 1:maxiter){     
  bb <- Tbb(par_Y, par_T, par_b, transf.par, dly, dlt, npar, GQ)
  wwN <- GetWeight(bb, par_Y, par_T, par_b, transf.par, dly, dlt, npar, GQ)
  Mstep <- M_step(bb, wwN, par_Y, par_T, par_b, transf.par,dly, dlt, npar)
  if(length(Mstep)==1){
	   err <- Mstep
	   break;
  }    
 	newpar_Y <- Mstep[1:nY]   	 
 	newpar_T <- Mstep[nY + 1:nTl]
 	newpar_b <- Mstep[(nY+nTl) + 1:nb]     
 	newlambda <- newpar_T[(nG+2):nTl]     
 	error <- c(newpar_Y-par_Y,newpar_T-par_T,newpar_b-par_b)    
 	n.conv <- rbind(n.conv, c(iter, max(abs(error)) ))  
 	if(max(abs(error))< epsilon){ 
		break 
 	}     
 	par_Y = newpar_Y;
 	par_T = newpar_T;                                   
 	par_b = newpar_b; 
}   # end of 'iter' #
 	
loglik=GetLogLik_phi_b1(newpar_Y,newpar_T,newpar_b, transf.par, dly,dlt,npar, GQ)  	
      BIC = -2*loglik + (nY+nTl+nb)*log(N) 
      AIC = -2*loglik + (nY+nTl+nb)*2  	

out=list(coef.lm.y=newpar_Y[1:nB], coef.nlm.y=newpar_Y[nB+1:n.basis], var.resid=newpar_Y[nY],
         raneff.vcomp=newpar_b,
         coef.fixed.surv=newpar_T[1:nG], coef.frailty.surv=newpar_T[nG+1:nP], lambdas= newpar_T[(nG+nP)+1:mt],
             loglik=loglik,AIC=AIC,BIC=BIC, 
         degree=degree, n.knots=n.knots, wt=dlt$wt, K=transf.par,
         vnf.y=varname.fixed.y, vnr.y=varname.random.y, vnf.t=varname.fixed.t,formula.y=formula(obj1))

if (CovEst==TRUE){
  COV <- covest(newpar_Y,newpar_T,newpar_b, transf.par, dly,dlt,npar, GQ)  
      errCOV= COV[[2]]; 
      COV=COV[[1]];
  var_Y <- COV[1:nY, 1:nY]
  var_T <- COV[nY+1:nTl, nY+1:nTl] 
  var_b <- COV[nY+nTl+1:nb, nY+nTl+1:nb]  
  out=list(coef.lm.y=newpar_Y[1:nB], coef.nlm.y=newpar_Y[nB+1:n.basis], var.resid=newpar_Y[nY],
           raneff.vcomp=newpar_b,
           coef.fixed.surv=newpar_T[1:nG], coef.frailty.surv=newpar_T[nG+1:nP], lambdas= newpar_T[(nG+nP)+1:mt],
               loglik=loglik,AIC=AIC,BIC=BIC, 
            degree=degree, n.knots=n.knots,  wt=dlt$wt, K=transf.par,
            vnf.y=varname.fixed.y, vnr.y=varname.random.y, vnf.t=varname.fixed.t, formula.y=formula(obj1),
            covy=var_Y,covb=var_b,covt=var_T)
}
class(out) <- 'jplm'
out
}