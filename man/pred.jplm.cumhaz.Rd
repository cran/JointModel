\name{pred.jplm.cumhaz}
\alias{pred.jplm.cumhaz}

\title{
Predict the baseline cumulative hazard function at any given time point
}

\description{
This function calculates a predicted baseline cumulative hazard function, evaluated at given time points.
}

\usage{
pred.jplm.cumhaz(object, at=NULL, CI=FALSE)
}

\arguments{
	\item{object}{a Joint Model fit object, i.e., the result of \code{jplm}.}	
	\item{at}{a vector of fixed time points to be evaluated.}
	\item{CI}{logical value; if \code{TRUE}, a 95\% pointwise confidence interval is returned. }
}

\value{
	If \code{CI=FALSE}, it returns a numeric vector of predicted cumulative hazard values at \code{at= }. If \code{CI=TRUE}, it returns a numeric vector of predicted cumulative hazard values, the standard error estimate of the predicted value, and its lower and upper 95\% pointwise confidence interval. 
}

\author{Sehee Kim}

\examples{

#	a simulated data set of longitudinal responses
attach(prostate)
#	a simulated data set of drop-out process (or, time-to-event)
attach(dropout)

#	joint fit of a partially linear model and a proportional odds model
#	with a subject-specific random intercept and random slope 
fit1 <- jplm(logPSA.postRT ~ logPSA.base + (1 + VisitTime|ID), 
                 nlm.par=prostate$VisitTime, data.y=prostate, 
                 Surv(DropTime, Status) ~ logPSA.base2, 
                 formula.frailty= ~ 1 + DropTime, 
                 id.vec=dropout$ID2, transf.par=1, data.surv=dropout)

# Evaluate at 20,...,80 percent of the maximum observed survival time 	
pts <- c(0.2, 0.4, 0.6, 0.8)*max(dropout$DropTime)  
pred.jplm.cumhaz(fit1, at=pts)
out <- pred.jplm.cumhaz(fit1, at=pts, CI=TRUE)
out$Value	
}

\references{
Kim, S., Zeng, D., Taylor, J.M.G. (2016) Joint partially linear model for longitudinal data with informative drop-outs. 
\emph{Under revision} \bold{0}, 000-000.
}

\seealso{
    \code{\link{jplm}}
}
