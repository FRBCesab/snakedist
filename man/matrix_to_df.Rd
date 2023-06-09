% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_to_df.R
\name{matrix_to_df}
\alias{matrix_to_df}
\title{Convert a distance matrix to a data.frame}
\usage{
matrix_to_df(x, all = FALSE)
}
\arguments{
\item{x}{a \code{matrix}.}

\item{all}{a \code{logical} value. If \code{TRUE} (default), keep cells with \code{NA}.}
}
\value{
A \code{data.frame} with the following three columns: \code{from} (the
first site), \code{to} (the second site) and \code{distance} (the value of the link
between two sites).
}
\description{
Converts a (distance) \code{matrix} to a \code{data.frame}.
}
\examples{
# Add an example ----
}
