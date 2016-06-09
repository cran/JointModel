
### Construct a sieve space of nonlinear trajectory
### KNOTS: use every 100/(Kn+1)th quantiles of observed Y measuring time;
basis <- 
function(nl.fixed, degree=3, n.knots=3,bdr){
  bsq= (1:n.knots)/(n.knots+1);      
  inner.knots= quantile(nl.fixed, prob=bsq)
  bs(nl.fixed,degree,inner.knots,intercept=TRUE,Boundary.knots=bdr)
}
