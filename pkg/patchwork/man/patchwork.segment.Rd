\name{patchwork.segment}
\alias{patchwork.segment}

\title{
patchwork segmentation function
}

\description{
Internal function of patchwork.plot().
}

\usage{
patchwork.segment(kbsegs,chroms,Alpha,SD)
}

\arguments{
  \item{kbsegs}{
  	kbsegs generated from patchwork.smoothing().
	}
  \item{chroms}{
	chroms is an object generated from pile.alleles.RData's alf\$achr attribute.
	}
	\item{Alpha}{
	Default is inherited from patchwork.plot() and is 0.0001.
	From DNAcopy (?segment):
	alpha: significance levels for the test to accept change-points.
	}
	\item{SD}{
	Default is inherited from patchwork.plot() and is 1.
	From DNAcopy (?segment):
	undo.SD: the number of SDs between means to keep a split if undo.splits="sdundo".
	}
}

\details{
Performs a segmentation using various functions of DNAcopy package and returns the segments
in a new object called "segs", saved in Segments.Rdata. \cr

See ?DNAcopy
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