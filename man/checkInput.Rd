% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CheckInputforPlotting.R
\name{checkInput}
\alias{checkInput}
\title{Return an error if input for DGObarplot and DGOnetplot is faulty.}
\usage{
checkInput(DGOResult)
}
\arguments{
\item{DGOResult}{a list of 2 enrichResults objects, DO and GO}
}
\value{
Returns an error if input is corrupted.
}
\description{
A helper function that checks the input of DGObarplot and DGOnetplot.
}
