\name{patchwork.applyref}
\alias{applyref}
\alias{patchwork.applyref}


\title{
Apply reference to query data
}

\description{
Applies reference data depending on the sequencing type of the data.
}

\usage{
patchwork.applyref(data,reference)
}

\arguments{
  \item{data}{
  	The data object generated from patchwork.GCNorm().
	}
  \item{reference}{
	Path to a reference file that you have either downloaded from us or created yourself with
	patchwork.createreference().
	}
}

\details{
This function applies a normal reference from Solexa (Illumina), SOLiD or
a reference that you have created for the data.
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
\code{\link{patchwork.createreference}}
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