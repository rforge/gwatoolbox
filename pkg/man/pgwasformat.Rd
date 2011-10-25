\name{pgwasformat}
\alias{pgwasformat}
\title{Formatting of GWAS result files.}
\description{Formats headers, orders columns, calculates inflation factors and applies genomic control in GWAS result files.}
\usage{pgwasformat(script, logfile, cpu)}
\arguments{
	\item{script}{
		Name of a textual input file with processing instructions. 
		The file should contain the names and locations of all GWAS data files to be processed 
		along with basic information from each individual study, and instructions for the quality check.
	}
	\item{logfile}{
		Name of a log file with processing output. 
	}
	\item{cpu}{
		Number of CPUs.
	}
}
\author{Daniel Taliun, Christian Fuchsberger, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{plot}
\keyword{package}
\examples{

}