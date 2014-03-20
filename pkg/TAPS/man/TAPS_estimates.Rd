\name{TAPS_estimates}
\alias{TAPS_estimates}
\alias{Taps_estimates}
\alias{TAPS.estimates}
\alias{Taps.estimates}


\title{
Tumor Aberration Prediction Suite - Estimate allele-specific copynumbers
}

\description{
TAPS_estimates() is intended to be used after TAPS_plot() and before TAPS_click().

TAPS_estimates() estimates the samples allele-specific copynumbers and loss of heterozygosity occurence.  
}


\usage{
TAPS_estimates(path=getwd())
}

\arguments{
	\item{path}{
	Default is getwd(). Specifying to a specific samples directory will run TAPS_estimates on that directory.
	Specifying a directory containing one or more subdirectories that are samples (and not any other subdirectories
	or TAPS_estimates will error when trying to run them!) will iteratively run TAPS_estimates on all samples.
	}
}

\details{
TAPS_estimates() automatically loads shortRegions.Rdata and attempts to estimate the allele-specific copynumbers. Note that it does not always do this correctly. As such it should be used primarily to simplify TAPS_plot output interpretation but should not be used INSTEAD OF TAPS_click().

For example:\cr
1. Run TAPS_plot().\cr
2. Run TAPS_estimates(). \cr
3. Run TAPS_click(), some of the estimates made by TAPS_estimates are good and we can accept them, some are bad and we choose to redo them.\cr
4. Run TAPS_call().\cr
}


\references{
http://www.biomedcentral.com/content/pdf/gb-2011-12-10-r108.pdf
}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Hanna Goransson Kultima, \email{hanna.goransson.kultima@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se} \cr
Bjorn Viklund, \email{bjorn.viklund@medsci.uu.se}
}
