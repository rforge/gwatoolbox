# Copyright © 2011 Daniel Taliun, Christian Fuchsberger and Cristian Pattaro. All rights reserved.
#
# This file is part of GWAtoolbox.
#
# GWAtoolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWAtoolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWAtoolbox.  If not, see <http://www.gnu.org/licenses/>.


dispersion_check <- function(summary_file, sample_sizes, plot) {

	summary_file <- as.character(summary_file)

       	filelist <- scan(summary_file, what="character", quiet=TRUE)
	f <- length(filelist)
	if (f <= 0) {
		stop("File '", summary_file, "' is empty.\n")
	}

	if (!is.null(sample_sizes)) {
		if (!is.vector(sample_sizes, mode="numeric")) {
			stop("Argument 'sample_sizes' must be a numeric vector of positive values.")
		}
		sample_sizes<-sample_sizes[!is.na(sample_sizes) & sample_sizes > 0]
		if (length(sample_sizes) != f) {
			stop("The length of vector 'sample_sizes' must be equal to the number of studies in 'summary_file'.")
		}
	}

	if (!is.logical(plot)) {
		stop("Argument 'plot' must be a logical.")
	}

        names <- c("study", "mean_se", "median_n");
	data <- data.frame(study=NA, mean_se=NA, median_n=NA)
	
	if (!is.null(sample_sizes)) {
		for (i in 1:f) {
			x <- read.table(filelist[i], sep=";", header=T, stringsAsFactors=F)
			data[i, 1] <- filelist[i]
			data[i, 2] <- x$STDERR[4]
			data[i, 3] <- sample_sizes[i]
		}
	}
	else {
		for (i in 1:f) {
			x <- read.table(filelist[i], sep=";", header=T, stringsAsFactors=F)
			data[i, 1] <- filelist[i]
			data[i, 2] <- x$STDERR[4]
			data[i, 3] <- x$N[8]
		}
	}

	if (plot) {
		par(mfrow = c(1, 1))
		plot(data$median_n, data$mean_se, xlab = "Sample Size", ylab = "Standard Error", col = "black", pch = 20, cex = 1)
	}

	return(data)
}