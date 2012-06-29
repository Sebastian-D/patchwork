\name{CG_KaChCN}
\alias{CG_KaChCN}
\alias{KaChCN}


\title{
Karyotype plotting for each chromosome post patchwork.CG.copynumbers().
}

\description{
Visualizes the calculated data of patchwork.CG.plot() + patchwork.CG.copynumbers() for each chromosome.
See details for a walkthrough of the plot.
}

\usage{
CG_KaChCN(chr,start,end,int,ai,Cn,mCn,mchr,mpos,mval,sval,schr,spos,name='',xlim=c(0,2.4),ylim=c(0,1),maxCn=8)
}

\arguments{
  	\item{chr}{
  	Chromosomes name of segm object.
	}
  	\item{start}{			
	Start position of segment.
	}
	\item{end}{					
  	End position of segment.	
	}
	\item{int}{					
  	Ratio of relative coverage from somaticCnvSegmentsNonDiplod file and segment length.
	}
	\item{ai}{					
  	Allelic imbalance.
	}
	\item{Cn}{
  	Copy number.
	}
	\item{mCn}{
  	mCn.
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
  	Name that you want output to have. Default is '' which will yield files named as '_KaChCN_chr#.png'.
  	When method is run by patchwork.CG.copynumbers() it will that functions name parameter.
	}
  	\item{xlim}{
  	Default is c(0,2.4). Limit of x axis.
	}
  	\item{ylim}{
  	Default is c(0,1). Limit of y axis.
	}
	\item{maxCn}{
	Default is 8. The highest copynumber of the sample to check for.
	}
}

\details{
Walkthrough of the plot:

\bold{TOP} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Total Intensity \cr

The chromosome plotted against the complete genome background. The separation between
clusters within the plot are due to the fluctuating intensity and allelic imbalance
and as such display the varying allele counts and copy numbers.
Longer/larger segments have bigger circles. Darker circles show more content as they are
ontop of eachother.

\bold{TOP MIDDLE} \cr
Vertical axis: Copynumber \cr
Horizontal axis: Chromosomal coordinate \cr

Displays the total and minor copynumbers for different segments of the chromosome in question.

\bold{LOWER MIDDLE} \cr
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
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.CG.copynumbers}}
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