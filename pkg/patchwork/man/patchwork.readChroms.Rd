\name{patchwork.readChroms}
\alias{patchwork.readChroms}
\alias{readChroms}


\title{
patchwork function to read chromosomes into R
}

\description{
Internal function of patchwork.plot().
}

\usage{
patchwork.readChroms(BamFile,chroms)
}

\arguments{
  \item{BamFile}{
  	Input Bam from patchwork.plot().
	}
  \item{chroms}{
	chroms is an object generated from pile.alleles.RData's alf\$achr attribute.
	}
}

\details{
Reads GC values regarding each chromosome and attaches it to the data. Also uses the python
script "Pysamloader.py" to generate a position file for each chromosome in turn.
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