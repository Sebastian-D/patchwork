\name{patchwork.copynumbers}
\alias{patchwork.copynumbers}
\alias{copynumbers}


\title{
Function plot allele-specific copy numbers of the genome
}

\description{
Input data for this function should be assessed from any of the karyotype_chroms()
plots. For example <yoursample>_karyotype.chr<1-22,X,Y>.png.

See argument details.

}

\usage{
patchwork.copynumbers(CNfile,cn2,delta,het,hom,maxCn=8,ceiling=1,forcedelta=F,male.sample=F,male2femref=F)
}

\arguments{
	\item{CNfile}{
	The name and path of your copynumbers file, generated from patchwork.plot(). Example Myfile_copynumbers.Rdata.
	}
  \item{cn2}{
	The approximate position of copy number 2,diploid, on total intensity / coverage axis.
	}
  \item{delta}{
	The difference in total intensity between consecutive copy numbers. For example 1 and 2 or 2 and 3.
	If copy number 2 has total intensity ~0.6 and copy number 3 har total intensity ~0.8 then delta would be 0.2.
	}
  \item{het}{
	Allelic imbalance ratio of heterozygous copy number 2.
	}
  \item{hom}{
	Allelic imbalance ratio of Loss-of-heterozygosity copy number 2.
	}
  \item{maxCn}{
	Highest copy number to calculate for. Default is 8.
	}
  \item{ceiling}{
	Default is 1.
	}
  \item{forcedelta}{
	Default is FALSE. If TRUE the delta value will not be subject to small adjustment changes.
	}
  \item{male.sample}{
	Default is FALSE. If it is a male sample put TRUE here and it will handle the XY chromosomes better.
	}
  \item{male2femref}{
	Default is FALSE. If TRUE the sample is male but the reference you used is female. This will correct for this.
	}
}

\details{
1. Run the patchwork.plot correctly. This will generate <yourfile>_copynumbers.Rdata.
2. View your plots and decide the correct parameters for the other arguments of this function.
4. Run the function.
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
\code{\link[patchwork]{patchwork.plot}} \cr
\code{\link[patchwork]{subfunctions_copynumbers}}
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