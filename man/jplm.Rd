\name{jplm}
\alias{jplm}

\title{
Joint Partially Linear Model for Longitudinal and Time-to-Event Data
}

\description{
This function fits jointly a partially linear model for normal longitudinal responses and a semiparametric transformation model for time-to-event data using a sieve maximum likelihood approach (Kim et al., 2016). 
}

\usage{
jplm(formula.lm.y, nlm.par=NULL, data.y, 
     formula.surv.fixed, formula.frailty, 
     id.vec=NULL, transf.par=0, data.surv, 
     degree=3, n.knots=3, CovEst=TRUE, maxiter=200, 
     epsilon=5e-04,\dots)
}

\arguments{
	\item{formula.lm.y}{a formula inheriting from class \code{lmer}. The right side of the formula describes both the linear fixed effects and random effects part of the longitudinal model. Random-effects terms are distinguished by one vertical bar (|) separating expressions for design matrices from grouping factors. Only one vertical var is allowed. }
	\item{nlm.par}{a vector of nonlinear effect covariate.}
	\item{data.y}{a data.frame containing the variables named in \code{formula.lm.y}. }
	\item{formula.surv.fixed}{a formula inheriting from class \code{coxph}. That is, the response must be a survival object as returned by the \code{Surv(,)} function, and the right side of formula must include only fixed effects covariates.}
	\item{formula.frailty}{a formula with the right side descring frailty term of the survival model. The left side of the formula must leave in blank, and the intercept term will be included unless specified by \code{ ~ -1+\dots} or \code{ ~ 0+\dots}. }
	\item{id.vec}{a vector containing subject ID corresponding to \code{formula.frailty}.}
	\item{transf.par}{a non-negative value of transformation parameter applied to the cumulative hazard function. \code{transf.par=0} will fit a proportional hazards model, while \code{transf.par=1} will fit a proportional odds model. Default is 0. }
	\item{data.surv}{a data.frame containing the variables named in the \code{formula.surv.fixed}.}
	\item{degree}{degree of the sieve polynomial. Default is 3 for cubic splines.}
	\item{n.knots}{the number of interior knot points for the B-spline approximation of the nonlinear effect. Default is 3.}
	\item{CovEst}{logical value; if \code{TRUE}, the covariate matrix of all the model parameters are estimated. }
	\item{maxiter}{the maximum number of EM iterations. Default is 200.}
	\item{epsilon}{tolerance value in the Newton-Raphson algorithm used to update the 
            parameters in the M-step. Default is 5e-04.}
    \item{\dots}{other arguments}
}

\value{
  	\item{coef.lm.y}{the vector of linear coefficients in the longitudinal model.}
  	\item{coef.nlm.y}{the vector of sieve coefficients corresponding to B-spline approximation of the nonlinear effect in the longitudinal model.}
  	\item{var.resid}{the variance estimate of longitudinal response residuals.}
  	\item{raneff.vcomp}{the vector of estimates of random effects variance component, corresponding to random intercept, random slope, and their correlation.}
  	\item{coef.fixed.surv}{the vector of fixed coefficients in the survival model.}
  	\item{coef.frailty.surv}{value of random (frailty) coefficient in the survival model.}
  	\item{lambdas}{the vector of jump sizes of the baseline cumulative hazard function, corresponding to the ordered observed event times. }
  	\item{loglik}{value of the log-likelihood with the final values of the coefficients.}
  	\item{AIC}{value of AIC with the final values of the coefficients.}
  	\item{BIC}{value of BIC with the final values of the coefficients.}
  	\item{degree}{degree used for the sieve polynomial.}
  	\item{n.knots}{number of interior knot points used for the B-spline approximation of the nonlinear effect.}
  	\item{K}{the transformation parameter.}
  	\item{covy}{the variance matrix of all the parameters in the longitudinal model, corresponding to the estimates of\code{coef.lm.y},\code{coef.nlm.y}, and \code{var.resid}.}
  	\item{covb}{the variance matrix of all the parameters in \code{raneff.vcomp}.} 		
  	\item{covt}{the variance matrix of all the parameters in the survival model, corresponding to the estimates of \code{coef.fixed.surv}, \code{coef.frailty.surv}, and \code{lambdas}.}
}

\details{
Function \code{jplm} fits joint models for longitudinal and survival data. Viewing time-to-drop-out as an event process, the same joint models also can be fitted to longitudinal data with informative drop-outs. For more detailed formulation of these models, refer to Kim et al. (2016).
   
For the longitudinal model specification, all linear effects part (including fixed and random)  should be described in  \code{formula.lm.y}, while a non-linear fixed effect component should be specified in the \code{nlm.par} argument. 
	
For the survival model specification, the fixed effects component should be described by the argument \code{formula.surv.fixed}, whereas the random effects (frailty) component should be described in \code{formula.frailty}. The current version assumes all random effects coefficients are the same. 
For example, the argument \code{formula.frailty= ~ 1+time} implies a regression model with \eqn{\phi(b_1 + b_2*time)}, where \eqn{b_1} and \eqn{b_2} represent the random intercept and slope terms and \eqn{\phi} is one dimensional coefficient corresponding to both \eqn{b_1} and \eqn{b_2}. 

A logarithmic transformation can be applied to the cumulative hazard function by varying the value of \code{transf.par}. Specifically, the transformation function takes the form of: \eqn{H(x) = 0} if \code{transf.par=0}; otherwise \eqn{H(x) = log(1 + K)/K} with \code{transf.par=K}. \code{transf.par=0} will fit a proportional hazards model, while \code{transf.par=1} will fit a proportional odds model.

For AIC- or BIC-based model selection, \code{CovEst=FALSE} is strongly recommended.
}

\author{Sehee Kim}

\examples{

# a simulated data set of longitudinal responses
attach(prostate)
# a simulated data set of time-to-event (e.g., drop-out process)
attach(dropout)

# joint fit of a partially linear model and a proportional hazards model 
# with a subject-specific random intercept and random slope
fit0 <- jplm(logPSA.postRT ~ logPSA.base + (1 + VisitTime|ID), 
                 nlm.par=prostate$VisitTime, data.y=prostate, 
                 Surv(DropTime, Status) ~ logPSA.base2, 
                 formula.frailty= ~ 1 + DropTime, 
                 id.vec=dropout$ID2, transf.par=0, data.surv=dropout)
summary(fit0)	
}

\references{
Kim, S., Zeng, D., Taylor, J. M. G. (2016) Joint partially linear model for longitudinal data with informative drop-outs. 
\emph{Under revision} \bold{0}, 000-000.
}

\seealso{
    \code{\link{pred.jplm.nonlinear}},
    \code{\link{pred.jplm.cumhaz}}
}
