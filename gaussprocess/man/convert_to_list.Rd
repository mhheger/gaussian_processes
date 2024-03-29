% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm_2.1.R
\name{convert_to_list}
\alias{convert_to_list}
\title{convert_to_list
A package interna, that provides the converting of the input data to the
needed list type}
\usage{
convert_to_list(x, n)
}
\arguments{
\item{x}{vector, matrix, data.frame, that should be converted into a list of
input vectors}

\item{n}{the dimension of each vector in the list that is returned}
}
\value{
a list of vectors of length n
}
\description{
convert_to_list
A package interna, that provides the converting of the input data to the
needed list type
}
\examples{
convert_to_list(matrix(1:6, nrow = 2), n= 2)
}
