\name{CG_KaCh}
\alias{CG_KaCh}
\alias{CG_KaryotypeChroms}

\title{
Karyotype plotting of each chromosome
}

\description{
Visualises the calculated data of patchwork.CG.plot() for each chromosome.
See details for a walkthrough of the plot.
}

\usage{
CG_KaCh(chr,start,end,int,ai,mchr,mpos,mval,sval,schr,spos,name,xlim,ylim)
}

\arguments{
  	\item{chr}{					
  	Chromosomes name of segm object, as generated from somaticCnvSegmentsNonDiplod file.
	}
  	\item{start}{				
	Start position of segment.
	}
	\item{end}{					
  	End position of segment.	
	}
	\item{int}{					
  	Ratio between relative coverage and segment length.
	}
	\item{ai}{					
  	Allelic imbalance as calculated from (1 - (lesser readcount / greater readcount)) per segment.
	}
	\item{mchr}{
	Chromosomes name of depcov (depthOfCoverage file).
	}
	\item{mpos}{
  	Start position from depcov.
	}
	\item{mval}{
  	Ratio between averagenormalizedcoverage and mean(averagenormalizedcoverage) from depthOfCoverage file.
	}
	\item{sval}{
	Inverse ratio between allele maximum and allele minimum per position.
	}
	\item{schr}{
	Chromosomes name of mastervar.
	}
	\item{spos}{
	Start position from mastervar.
	}
	\item{name}{
  	Name that you want output to have. Default is ' ' which will yield files named as '_KaCh_chr#.png'.
  	When method is run by patchwork.CG.plot() it will inherit the name parameter thereof.
	}
  	\item{xlim}{
  	Default is c(0,2.4). Limit of x axis for whole genome plot.
	}
  	\item{ylim}{
  	Default is c(0,1). Limit of y axis for whole genome plot.
	}
}

\details{
Walkthrough of the plot: \cr

\bold{TOP} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Total Intensity \cr

The chromosome plotted against the complete genome background. The separation between
clusters within the plot are due to the fluctuating intensity and allelic imbalance
and as such display the varying allele counts and copy numbers.
Longer/larger segments have bigger circles. The varying coloration of the segments correspond to
the amount of segments, since greater color scale variation is needed for a greater number of segments,
and where they are located on the chromosome. They coincide and should be compared with the middle plot.
Darker grey circles show more content as they are ontop of eachother.

\bold{MIDDLE} \cr
Vertical axis: Total Intensity \cr
Horizontal axis: Chromosomal coordinate \cr

The chromosome in questions total intensity plotted against the position on the chromosome.

\bold{BOTTOM} \cr
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
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.CG.plot}} \cr
\code{\link{CG_karyotype}}
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%Run examples of karyotype stuff. or not, maybe just point back to CNA.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line