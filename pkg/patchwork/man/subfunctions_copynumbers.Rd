\name{subfunctions_copynumbers}
\alias{subfunctions_copynumbers}
\alias{patchwork.subfunctions_copynumbers}
\alias{weightedMedian}
\alias{weightedMean}
\alias{is.autosome}
\alias{deChrom_ucsc}
\alias{chrom_ucsc}
\alias{assignAI}

\title{
A number of subfunctions used within patchwork.copynumbers().
}

\description{
A number of subfunctions used within patchwork.copynumbers().
These functions are identical in patchworkCG package, so should that package be installed
before patchwork or vice versa you will get a non-consequential message saying:

The following object(s) are masked from 'package:patchworkCG':

    chrom_ucsc, deChrom_ucsc, is.autosome, weightedMean, weightedMedian
    
    
}

\usage{
weightedMedian(data,weights)
weightedMean(data,weights)
is.autosome(vector)
deChrom_ucsc(data)
chrom_ucsc(data)
assignAI(snp_chr,snp_pos,snp_min,snp_max,seg_chr,seg_start,seg_end)
}

\arguments{
  \item{data}{
  	Object to be handled by the internal function.
	}
  \item{weights}{
	Weights to be applied.
	}
  \item{vector}{
	The segs objects chr column.
	}
\item{snp_chr}{
	Chromosome.
	}
\item{snp_pos}{
	Position.
	}
\item{snp_min}{
	Lesser allele amount in single position.
	}
\item{snp_max}{
	Greater allele amount in single position.
	}
\item{seg_chr}{
	Segments chromosome.
	}
\item{seg_start}{
	Start position of segment.
	}
\item{seg_end}{
	End position of segment.
	}
}

%\details{
%Details
%}

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
\code{\link{patchwork.copynumbers}}
}

%\examples{
%\dontrun{}
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%Run examples of karyotype stuff. or not, maybe just point back to patchwork.plot.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line