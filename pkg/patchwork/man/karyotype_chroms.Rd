\name{karyotype_chroms}
\alias{karyotype_chroms}

\title{
Karyotype plotting of each chromosome
}

\description{
Visualises the calculated data of patchwork.plot() for each chromosome.
See details for a walkthrough of the plot.
}

\usage{
karyotype_chroms(chr,start,end,int,ai,mchr,mpos,mval,schr,spos,sval,name,xlim,ylim)
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
  	Median.
	}
	\item{ai}{					%segs$ai
  	Allelic imbalance.
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
	\item{name}{
  	Name of sample. patchwork.plot sends it as name of the bamfile.
	}
  	\item{xlim}{
  	Default is c(0,2.4). Limit of x axis.
	}
  	\item{ylim}{
  	Default is c(0.3,1). Limit of y axis.
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
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.plot}} \cr
\code{\link{karyotype}}
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