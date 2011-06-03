\name{kusk_check}
\alias{kusk_check}
\title{Comparison Of Skewness and Kurtosis of Effect between GWAS studies.}
\description{Plots Skewness vs Kurtosis of Effect in GWAS studies.}
\usage{kusk_check(script, worst = c(50, 75, 90, 99, 100), plot = TRUE)}
\arguments{
	\item{script}{
		Name of a textual input file with processing instructions. 
		The file should contain the names and locations of all GWAS data files to be processed 
		along with basic information from each individual study, and instructions for the quality check.
	}
	\item{worst}{Percentage of worst p-values.}
	\item{plot}{Produce plot.}
}
\note{Under development.}
\author{Daniel Taliun, Christian Fuchsberger, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}