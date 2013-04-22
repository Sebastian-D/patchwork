\name{TAPS_region}
\alias{TAPS_region}
\alias{Taps_region}
\alias{TAPS.region}
\alias{Taps.region}


\title{
Tumor Abberation Prediction Suite - Region plotter
}

\description{
TAPS_region allows you to look closer at a region of the genome that you specify. It also displays
the known genes of the region.
}

\usage{
TAPS_region(directory=NULL,chr,region,hg18=F)
}

\arguments{
	\item{directory}{
	Default is getwd().
	}
	\item{chr}{
	Supply the chromosome you wish to look at. ex: chr=1
	}
	\item{region}{
	The region of the chromosome you are interested in. ex: region=10000000:13000000
	}
	\item{hg18}{
	Default is FALSE, meaning the hg19 known gene list will be used. If TRUE the hg18 known gene list is used.
	}
}

\details{
TAPS_region takes the directory of the sample that you want to look at (defaults to current working directory
if no input is given), the chromosome and the actual chromosomal coordinates as input to give you an in-depth
view.
Unlike any of the other plots this view also shows known genes from a gene list of your choice, hg18 or hg19.
}

\references{
http://www.biomedcentral.com/content/pdf/gb-2011-12-10-r108.pdf
}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Hanna Goransson Kultima, \email{hanna.goransson.kultima@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}

