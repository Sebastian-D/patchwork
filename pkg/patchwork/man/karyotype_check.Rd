\name{karyotype_check}
\alias{karyotype_check}
\alias{check}


\title{
Karyotype check function
}

\description{
Plots the whole genome coverage vs allelic imbalance with the approximated
areas copynumbers and allele constitution. Karyotype_check is called from 
patchwork.copynumbers function.
}

\usage{
karyotype_check(chr,start,end,int,ai,Cn,mCn,t,name='',xlim=c(-1.02,1.02),ylim=0:1)
}

\arguments{
  	\item{chr}{					%segs$chr
  	Chromosome.
	}
  	\item{start}{				%segs$start
	Start position of segment.
	}
	\item{end}{					%segs$end
  	End position of segment.	
	}
	\item{int}{					%segs$mean
  	Mean.
	}
	\item{ai}{					%segs$ai
  	Allelic imbalance.
	}
	\item{Cn}{
  	Copy number.
	}
	\item{mCn}{
  	mCn.
	}
	\item{t}{
  	A list between int and ai.
	}
	\item{name}{
  	Default is BamFile argument from patchwork.plot is used as argument.
	}
  	\item{xlim}{
  	Default is c(-1.02,1.02). Limit of x axis.
	}
  	\item{ylim}{
  	Default is 0:1. Limit of y axis.
	}
}


\details{
Vertical axis: Allelic Imbalance. \\
Horizontal axis: Relative coverage. \\ %Ask Marcus about this.

The naming scheme is Copynumber-m-LesserAlleleDistribution
So for example 2m0 means copynumber = 2, both alleles are the same
whereas 2m1 means copynumber = 2, 1 allele each.

Another example:
4m0, copynumber = 4, All allels are the same. (Loss of heterozygosity).
4m1, copynumber = 4, 3 alleles are the same, one is different.
4m2, copynumber = 4, 2 alleles each.

The total number of alleles present are always the copynumber.

Run within patchwork.copynumbers function. copynumbers.Rdata must be in working directory for that function
to run, and thus this one aswell as most of its input is taken from segs.
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
\code{\link{patchwork.copynumbers}}
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