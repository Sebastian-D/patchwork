\name{karyotype_chromsCN}
\alias{karyotype_chromsCN}
\alias{chromsCN}


\title{
Karyotype plotting for each chromosome post patchwork.copynumbers().
}

\description{
Visualises the calculated data of patchwork.plot() + patchwork.copynumbers() for each chromosome.
See details for a walkthrough of the plot.
}

\usage{
karyotype_chromsCN(chr,start,end,int,ai,Cn,mCn,mchr,mpos,mval,schr,spos,sval,t,name='',parameters,xlim=c(0,2.5),ylim=0:1,maxCn=8)
}

\arguments{
  	\item{chr}{					%segs$chr
  	Chromosomes name of segs object. Final version of segs created from patchwork.Medians_n_AI().
	}
  	\item{start}{				%segs$start
	Start position of segment.
	}
	\item{end}{					%segs$end
  	End position of segment.	
	}
	\item{int}{					%segs$mean
  	Normalized coverage.
	}
	\item{ai}{					%segs$ai
  	Allelic imbalance.
	}
	\item{Cn}{
  	Copy number.
	}
	\item{mCn}{
  	mCn.
	}
	\item{mchr}{
	Chromosomes name of kbsegs object. Final version of kbsegs created from patchwork.smoothing().
	}
	\item{mpos}{
  	Middle position of kbsegs segments.
	}
	\item{mval}{
  	Ratio between kbsegs coverage and refcoverage.
	}
	\item{schr}{
  	Chromosome names of alf object. Final version of alf created from patchwork.alleledata().
	}
	\item{spos}{
  	Positions in alf.
	}
	\item{sval}{
	Inverse ratio between allele maximum and allele minimum.
	}
	\item{t}{
	A list of approximate copynumber and allelic imbalance positions.
	}
	\item{name}{
  	Name of sample. patchwork.plot sends it as name of the bamfile.
	}
	\item{parameters}{
	The parameters sent to patchwork.copynumbers() to be plotted.
	}
  	\item{xlim}{
  	Default is c(0,2.5). Limit of x axis.
	}
  	\item{ylim}{
  	Default is 0:1. Limit of y axis.
	}
	\item{maxCn}{
	Default is 8. The highest copynumber of the sample to check for.
	}
}

\details{
Walkthrough of the plot:

\bold{LEFT} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Total Intensity \cr

The chromosome plotted against the complete genome background. The separation between
clusters within the plot are due to the fluctuating intensity and allelic imbalance
and as such display the varying allele counts and copy numbers.
Longer/larger segments have bigger circles. Darker circles show more content as they are
ontop of eachother.

\bold{TOP RIGHT} \cr
Vertical axis: Copynumber \cr
Horizontal axis: Chromosomal coordinate \cr

Displays the total and minor copynumbers for different segments of the chromosome in question.

\bold{UPPER MIDDLE RIGHT} \cr
Vertical axis: Total Intensity \cr
Horizontal axis: Chromosomal coordinate \cr

The chromosome in questions total intensity plotted against the position on the chromosome.

\bold{LOWER MIDDLE RIGHT} \cr

Cytoband information.

\bold{BOTTOM RIGHT} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Chromosomal coordinate \cr

The chromosome in questions allelic imbalance plotted against the position on the chromosome.
}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.copynumbers}}
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%Run examples of karyotype stuff. or not, maybe just point back to patchwork.plot.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line