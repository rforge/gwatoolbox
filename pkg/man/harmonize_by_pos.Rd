\name{harmonize_by_pos}
\alias{harmonize_by_pos}
\title{The harmonization of SNP identifiers and alleles.}
\description{Function harmonizes SNP identifiers and alleles in the input file according to the provided VCF reference.}
\usage{harmonize_by_pos(input, output, vcf, id = "markername", chromosome = "chr", position="position", alleles = c("first_allele", "second_allele"), flip = TRUE, sep = "\t", drop = FALSE, gzip = TRUE)}
\arguments{
	\item{input}{
		Name of an input file.
	}
	\item{output}{
		Name of an output file.
	}
	\item{vcf}{
		Name of the VCF reference.
	}
	\item{id}{
		Name of the column containing SNP identifiers in the specified input file.
	}
	\item{chromosome}{
		Name of the column containing chromosome names in the specified input file.
	}
	\item{position}{
		Name of the column containing chromosomal positions in the specified input file.
	}
	\item{alleles}{
		Vector with column names containing first and second alleles in the specified input file.
	}
	\item{flip}{
		If TRUE, automatically flips alleles in order to match variant strand in VCF reference.
		Otherwise, the variant is considered as unresolved.
	}
	\item{sep}{
		The field separator character in the specified input file.
	}
	\item{drop}{
		If TRUE, then all unresolved variants are dropped from the output file.
	}
	\item{gzip}{
		TRUE if output file is in gzip format.
	}
}
\author{Daniel Taliun, Christian Fuchsberger, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}