\docType{data}
\name{prostate}
\alias{prostate}
\title{Prostate Cancer Data: Part 1 - a simulated example of longitudinal responses}
\format{A data frame with 697 observations (n = 100 patients) on the following 4 variables.
  \describe{
    \item{\code{ID}}{a numeric vector of patient  ID.}
    \item{\code{logPSA.postRT}}{a numeric vector containing Prostate-specific Antigen (PSA) levels after radiation therapy, i.e., log(PSA(t) + 0.1) observed at time t.}
    \item{\code{VisitTime}}{a numeric vector of visiting time.}
    \item{\code{logPSA.base}}{a numeric vector of log(baseline PSA + 0.1).}
  }
}
\description{
  A simulated data set used in Kim et al. (2016) to illustrate the JPLM method. This data set was generated under settings mimicking the prostate cancer study.  For detailed data generation settings, see Kim et al. (2016).	
}
\source{ 
  The prostate cancer data have been previously studied by Proust-Lima et al. (2008) and Taylor et al. (2013), among others.
}
\references{
Kim, S., Zeng, D., Taylor, J. M. G. (2016) Joint partially linear model for longitudinal data with informative drop-outs. 
\emph{Under revision} \bold{0}, 000-000.

Proust-Lima, C., Taylor, J. M. G., Williams, S. G., Ankerst, D. P., Liu, N., Kestin, L. L., Bae, K., and Sandler, H. M. (2008) Determinants of change in prostate-specific antigen over time and its association with recurrence after external beam radiation therapy for prostate cancer in five large cohorts. \emph{International Journal of Radiation Oncology Biology Physics} \bold{72}, 782-791.

Taylor, J. M. G., Part, Y., Ankerst, D. P., Proust-Lima, C., Williams, S., Kestin, L., Bae, K., Pickles, T., and Sandler, H. (2013) Real-time individual predictions of prostate cancer recurrence using joint models. \emph{Biometrics} \bold{69}, 206-213.
}
\keyword{datasets}