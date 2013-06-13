\name{harmonize}
\alias{harmonize}
\title{The harmonization of SNP identifiers and alleles.}
\description{Function harmonizes SNP identifiers and alleles in the input file according to the provided VCF reference.}
\usage{harmonize(input, output, vcf, chromosome = "chr", id = "markername", alleles = c("ref_allele", "non_ref_allele"), sep = "\t", vcf_alleles = TRUE, drop = FALSE, gzip = TRUE)}
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
	\item{chromosome}{
		Name of the column containing chromosome names in the specified input file.
	}
	\item{id}{
		Name of the column containing SNP identifiers in the specified input file.
	}
	\item{alleles}{
		Vector with column names containing reference and non-reference alleles in the specified input file.
	}
	\item{sep}{
		The field separator character in the specified input file.
	}
	\item{vcf_alleles}{
		When applicable, alleles in the input file are recoded to R(reference)/I(insertion)/D(deletion) format.
		If vcf_alleles is set to TRUE, then alleles are recoded based on the VCF reference.
		Otherwise, alleles are recoded using only information available in the input file.
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