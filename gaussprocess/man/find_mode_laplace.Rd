% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification_laplace.R
\name{find_mode_laplace}
\alias{find_mode_laplace}
\title{finding value of mode function using binary laplace-approximation}
\usage{
find_mode_laplace(K, y, likelihood_fun)
}
\arguments{
\item{K}{covariance-variance-matrix of dimension (n,n)}

\item{y}{output-vector of length n, its values are just allowed
to be -1 or 1}

\item{likelihood_fun}{name of the response function, that should be used
either 'probit' or 'logit'}
}
\value{
mode: value of the mode function
         log_marginal_likelihood: logarithmic marginal likelihood of the
         result
}
\description{
finding value of mode function using binary laplace-approximation
}
\examples{
find_mode_laplace(  K = matrix( 1:4, nrow= 2),
                    y = c(1,-1),
                    likelihood_fun = "probit"
                    )
}
