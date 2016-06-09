\name{pred.jplm.nonlinear}
\alias{pred.jplm.nonlinear}

\title{
Predict a smoothed nonlinear effect on the longitudinal response
}

\description{
This function calculates a predicted nonlinear effect function evaluated at given points.
}

\usage{
pred.jplm.nonlinear(object, nlm.par, at=NULL, CI=FALSE)
}

\arguments{
	\item{object}{a Joint Model fit object, i.e., the result of \code{jplm}.}	
	\item{nlm.par}{a vector of nonlinear effect covariate, as specified in \code{nlm.par= } of \code{jplm}.}
	\item{at}{a vector of fixed points to be evaluated.}
	\item{CI}{logical value; if \code{TRUE}, a 95\% pointwise confidence interval is returned. }
}

\value{
	If \code{CI=FALSE}, it returns a numeric vector of predicted nonlinear effect at \code{at= }, the standard error estimate of the predicted value, and test result based on the asymptotic normality. If \code{CI=TRUE}, it returns a numeric vector of predicted nonlinear effect, the standard error estimate of the predicted value, and its lower and upper 95\% pointwise confidence interval. 
}

\author{Sehee Kim}

\examples{

#	a simulated data set of longitudinal responses
attach(prostate)
#	a simulated data set of time-to-event (e.g., drop-out process)
attach(dropout)

#	joint fit of a partially linear model and a proportional odds model
#	with a subject-specific random intercept and random slope 
fit1 <- jplm(logPSA.postRT ~ logPSA.base + (1 + VisitTime|ID), 
                 nlm.par=prostate$VisitTime, data.y=prostate, 
                 Surv(DropTime, Status) ~ logPSA.base2, 
                 formula.frailty= ~ 1 + DropTime, 
                 id.vec=dropout$ID2, transf.par=1, data.surv=dropout)

# Evaluate at 20,...,80 percent of the maximum measurement time 	
pts <- c(0.2, 0.4, 0.6, 0.8)*max(prostate$VisitTime)         
pred.jplm.nonlinear(fit1, prostate$VisitTime, at=pts)
out <- pred.jplm.nonlinear(fit1, prostate$VisitTime, at=pts, CI=TRUE)
out$Value	
}

\references{
Kim, S., Zeng, D., Taylor, J.M.G. (2016) Joint partially linear model for longitudinal data with informative drop-outs. 
\emph{Under revision} \bold{0}, 000-000.
}

\seealso{
    \code{\link{jplm}}
}
