\name{patchwork.Medians_n_AI}
\alias{patchwork.Medians_n_AI}
\alias{Medians_n_AI}


\title{
patchwork function which computes medians and allelec imbalance of the segment.
}

\description{
Computes medians and Allelic imbalance for segments of the input genome.
}

\usage{
patchwork.Medians_n_AI(segs,kbsegs,alf)
}

\arguments{
  \item{segs}{
  	Generated from patchwork.segment().
	}
  \item{kbsegs}{
	Generated from patchwork.smoothing().
	}
  \item{alf}{
	Generated from patchwork.alleledata().
	}
}

\details{
See homepage at http://patchwork.r-forge.r-project.org/ for more information.
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
\code{\link{patchwork.plot}}
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