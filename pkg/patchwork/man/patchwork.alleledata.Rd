\name{patchwork.alleledata}
\alias{alleledata}
\alias{patchwork.alleledata}


\title{
patchwork alleledata function
}

\description{
Loads the data object commonSnps132.RData and uses the perl script, included in package,
pile2alleles.pl to generate pile.alleles in your working directory. This will be read into R
as the alf object and then the file pile.alleles removed.
}

\usage{
patchwork.alleledata(Pileup, normalalf=NULL)
}

\arguments{
  \item{Pileup}{
  	Pileup file from patchwork.plot() input.
	}
\item{normalalf}{
	normalalf is generated and used if you have supplied a normal.pileup in its argument
	at patchwork.plot.
}
}

\details{
Uses SNP data to calculate allele frequencies of the sample. Currently uses SNP132.

patchwork.alleledata is a subfunction to patchwork.plot.
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