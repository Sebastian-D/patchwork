\name{patchwork.region}
\alias{Patchwork.Region}
\alias{region}
\alias{patchwork.region}


\title{
Detailed view of a selected regions of the genome.
}

\description{
Visualises a selected regions on a selected chromosome for closer inspection. Also shows known genes for selected region.
}

\usage{
patchwork.region(CNfile=NULL,chr,region,hg18=F)
}

\arguments{
  	\item{CNfile}{
  	The <sample>_copnumbers.Rdata file generated from patchwork.plot(). If kept at default, NULL, the current working directory will be searched for it.
	}
  	\item{chr}{	
	Which chromosome you wish to view a region of. ex chr="1"
	}
	\item{region}{
  	The region in chromosomal coordinates you wish to view. ex region=1000000:6000000
	}
	\item{hg18}{
  	Gene list for hg19 aligned genome is default. Set hg18=T to use hg18 gene list.
	}
}

\details{
Walkthrough of the plot:

\bold{TOP LEFT} \cr
Vertical axis: Allelic imbalance \cr
Horizontal axis: Normalized coverage \cr

The selected regions segments in color, plotted agains the whole genome background in grey.

\bold{TOP RIGHT - TOP} \cr
Vertical axis: Normalized coverage \cr
Horizontal axis: Position (Mb) \cr

Normalized coverage of chromosome with selected region highlighted in a transparent red bar.

\bold{TOP RIGHT - MIDDLE} \cr
Horizontal axis: Position (Mb) \cr

Cytoband information for chromosome with selected region highlighted in a transparent red bar.

\bold{TOP RIGHT - BOTTOM} \cr
Vertical axis: Allelic imbalance \cr
Horizontal axis: Position (Mb) \cr

Allelic imbalance of chromosome with selected region highlighted in a transparent red bar.

\bold{BOTTOM - TOP} \cr
Vertical axis: Normalized coverage \cr
Horizontal axis: Position (Mb) \cr

Normalized coverage of the selected region.

\bold{BOTTOM - MIDDLE} \cr
Horizontal axis: Position (Mb) \cr

Known gene positions which are in the selected region.

\bold{BOTTOM - BOTTOM} \cr
Vertical axis: Allelic imbalance \cr
Horizontal axis: Position (Mb) \cr

Allelic imbalance of the selected region.

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

% \seealso{
% \code{\link{patchwork.copynumbers}}
% }

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