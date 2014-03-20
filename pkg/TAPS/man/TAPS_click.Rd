\name{TAPS_click}
\alias{TAPS_click}
\alias{Taps_click}
\alias{TAPS.click}
\alias{Taps.click}


\title{
Tumor Aberration Prediction Suite - "GUI" simplified plot interpretation 
}

\description{
TAPS_click is intended to be used after \code{\link{TAPS_plot}} and TAPS_estimates and before TAPS_call.

TAPS_click is a help function designed for faster interpretation of TAPS_plot results. It iterates through all available samples, produces an interactive GUI and automatically forwards the user interpretation to SampleData.csv.
}


\usage{
TAPS_click(path = getwd())
}

\arguments{
	%\item{samples}{
	%Default is 'all'. Alternatively, supply a (boolean, integer or sample name) vector of included samples.
	%}
	\item{path}{
	Default is getwd(). Specifying to a specific samples directory will run TAPS_click on that directory.
	Specifying a directory containing one or more subdirectories that are samples (and not any other subdirectories
	or TAPS_click will error when trying to run them!) will iteratively run TAPS_click on all samples.
	}
}

\details{
TAPS_click takes the directory of the samples that you want to look at (defaults to the current working directory if no input is given) and iterates through all subdirectories. Each subdirectory is assumed to be a sample. For each sample it produces a GUI for the user showing the TAPS_plot generated plot and waits for the user to provide the interpretation. Upon completion it will send the user interpretation input to SamplesData.csv which is the file that TAPS_call uses for it's allele-specific copynumber analysis.
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
