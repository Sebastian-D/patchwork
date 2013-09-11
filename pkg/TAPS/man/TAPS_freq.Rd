\name{TAPS_freq}
\alias{TAPS_freq}
\alias{TAPS.freq}
\alias{Taps.freq}
\alias{Taps_freq}


\title{
Tumor Abberation Prediction Suite
}

\description{
TAPS_freq generates a summary of aberration frequencies in several samples (on which TAPS_call has completed).
}

\usage{
TAPS_freq(samples='all', outdir='frequencies', hg19=T)
}

\arguments{
	\item{samples}{
	Default is 'all'. Alternatively, supply a (boolean, integer or sample name) vector of included samples.
	}
	\item{outdir}{
	All results are written to this directory.
	}
	\item{hg19}{
	The reference genome 'hg19' will be used for genes. If FALSE, reverts to hg18.
	}
}

\details{
This function summarizes aberration frequencies over multiple samples. 
Run from the directory containing sample folders.
Make sure TAPS_call() completed successfully on all (included) samples.
}

\references{
http://www.biomedcentral.com/content/pdf/gb-2011-12-10-r108.pdf
}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Hanna Goransson Kultima, \email{hanna.goransson.kultima@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
