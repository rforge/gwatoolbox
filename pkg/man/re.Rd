\name{re}
\alias{re}
\title{Random Effects Heterogeneity Analysis.}
\description{Outlier cohorts detection in meta-analysis based on random effects.}
\usage{re(script)}
\arguments{
	\item{script}{
		Name of a textual input file with processing instructions. 
		The file should contain the names and locations of all GWAS data files to be processed 
		along with basic information from each individual study.
	}
}
\section{Specifying The Input Data Files}{
	The names of the GWAS data files are specified in the input script with the command \bold{PROCESS} (one line per file). 
	A different directory path can be specified for each file.
	
	Example:
	\tabular{l}{
		PROCESS   input_file_1.txt\cr
		PROCESS   /dir_1/dir_2/input_file_2.csv	
	}
	The annotation is applied first to \file{input_file_1.txt} and then to \file{input_file_2.csv}.
}
\section{Specifying Column Names in Input Data Files}{
	In the table below, the complete list of the default column names for the GWAS data file is reported. 
	These names identify uniquely the items in the GWAS data file.
	\tabular{ll}{
		Default column name(s)\tab Description\cr
		MARKER\tab Marker name\cr
		CHR\tab Chromosome number or name\cr
		POSITION\tab Marker position
	}
	Given that different names can be provided for each GWAS data file, \code{annotate()} allows to redefine the default values for 
	every input file in the input script. 
	The redefinition command consists of the default column name followed by the new column name.
	When the map file is specified using command \bold{MAP_FILE}, then \bold{CHR} and \bold{POSITION} columns in the GWAS data file are not required.
	
	Example:

	Let's assume to have two input files, \file{input_file_1.txt} and \file{input_file_2.csv}. 
	In the \file{input_file_1.txt}, the column names for marker name, chromosome name and position are \emph{SNPID}, \emph{CHR} and \emph{POS}, respectively. 
	In the \file{input_file_2.csv}, the column names for marker name is the same as in \file{input_file_1.txt}, but the column names 
	for the chromosome and position are \emph{chromosome} and \emph{position}, respectively. 
	The correct column redefinition is as follows:
	\tabular{l}{
		MARKER     SNPID\cr
		POSITION   POS\cr
		PROCESS    input_file_1.txt\cr
		CHR        chromosome\cr
		POSITION   position\cr
		PROCESS    input_file_2.csv
	}
	There are no need to define the \bold{CHR} field for the \file{input_file_1.txt}, since it matches the default name.
}
\section{Field Separator in Input Data Files}{
	The field (column) separator can be different for each GWAS data file.
	\code{annotate()} automatically detects the original separator field for each input file \emph{based on the first 10 rows}. 
	However, the user has the possibility to specify the original separator manually for each individual file using the command \bold{SEPARATOR}. 
	The supported arguments and related separators are listed below:
	\tabular{ll}{
		Argument\tab Separator\cr
		COMMA\tab \emph{comma}\cr
		TAB\tab \emph{tabulation}\cr
		WHITESPACE\tab \emph{whitespace}\cr
		SEMICOLON\tab \emph{semicolon}
	}

	Example:
	\tabular{l}{
		PROCESS      input_file_1.txt\cr
		SEPARATOR    COMMA\cr              
		PROCESS      input_file_2.csv\cr
		PROCESS      input_file_3.txt 
	}
	For the input file \file{input_file_1.txt} the field separator is determined automatically by the program but, 
	for the input files \file{input_file_2.csv} and \file{input_file_3.txt} the separator is manually set to comma by the user.
}
\section{Case Sensitivity}{
	By default the \code{annotate()} assumes that column names in the all specified files are case insensitive. 
	For example, the column names \emph{CHR}, \emph{Chr}, and \emph{chr} are all perfectly equivalent.
	This behaviour can be modified for every input file in the input script using the command \bold{CASESENSITIVE}, 
	that controls case sensitivity for the column names, as specified below:
	\tabular{ll}{
		Argument\tab Description\cr
		0\tab Column names in the input file\cr
		\tab are case insensitive (default)\cr
		1\tab Column names in the input file\cr
		\tab are case sensitive
	}

	Example:
	\tabular{l}{
		CASESENSITIVE   1\cr
		MAP_FILE        map_file_1.txt\cr
		REGIONS_FILE    genes_file_1.txt\cr
		PROCESS         input_file_1.txt\cr
		CASESENSITIVE   0\cr
		MAP_FILE        map_file_2.csv\cr
		REGIONS_FILE    genes_file_2.csv\cr
		PROCESS         input_file_2.csv 
	}	
}
\section{Output File Name}{
	The output file names are created by adding a prefix to the input file names.	
	The prefix is specified with the command \bold{PREFIX}. 
	
	Example:
	\tabular{l}{
		REGIONS_FILE genes_file_1.txt\cr
		PREFIX       annotated_\cr
		PROCESS      input_file_1.txt\cr
		PROCESS      input_file_2.csv\cr
		PREFIX       annot_\cr
		PROCESS      input_file_3.tab\cr
	}
	All the output files corresponding to the input files \file{input_file_1.txt} and \file{input_file_2.csv} will be prefixed with \emph{annotated_}; 
	the output files corresponding to the input file \file{input_file_3.tab} will be prefixed with \emph{annot_}.
}
\author{Daniel Taliun, Christian Fuchsberger, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}
\examples{}