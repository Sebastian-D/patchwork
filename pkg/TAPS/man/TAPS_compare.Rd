\name{TAPS_compare}
\alias{TAPS_compare}
\alias{TAPS_compare}
\alias{TAPS_compare}
\alias{TAPS_compare}


\title{
Tumor Abberation Prediction Suite
}

\description{
TAPS_compare generates a comparison of aberration frequencies in two groups of samples (on which TAPS_call has completed).
}

\usage{
TAPS_compare(grp1, grp2, name1='1', name2='2', outdir='frequencies_comp', hg19=T)
}

\arguments{
	\item{grp1}{
	Supply a (boolean, integer or sample name) vector of samples in the first group.
	}
	\item{grp2}{
	Supply a (boolean, integer or sample name) vector of samples in the second group.
	}
	\item{name1}{
	Supply a name for the first group. (Default='1')
	}
	\item{name2}{
	Supply a name for the second group. (Default='2')
	}
	\item{outdir}{
	All results are written to this directory.
	}
	\item{hg19}{
	The reference genome 'hg19' will be used for genes. If FALSE, reverts to hg18.
	}
}

\details{
This function compares aberration frequencies in two groups of samples. 
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
