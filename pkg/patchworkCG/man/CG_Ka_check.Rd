\name{CG_Ka_check}
\alias{CG_Ka_check}
\alias{check}


\title{
Karyotype check function
}

\description{
Plots the whole genome coverage vs allelic imbalance with the approximated
areas copynumbers and allele constitution. CG_Ka_check is called from the patchwork.CG.copynumbers() function.
It would be hard to run this function outside of patchwork.CG.copynumbers().
}

\usage{
CG_Ka_check(chr,start,end,int,ai,Cn,mCn,t,name='',xlim=c(0,2.4),ylim=c(0,1))
}

\arguments{
  	\item{chr}{	
  	Chromosome.
	}
  	\item{start}{
	Start position of segment.
	}
	\item{end}{	
  	End position of segment.	
	}
	\item{int}{	
  	Ratio between relative coverage and segment length.
	}
	\item{ai}{
  	Allelic imbalance as calculated from (1 - (lesser readcount / greater readcount)) per segment.
	}
	\item{Cn}{
  	Copy number object created in patchwork.CG.copynumbers().
	}
	\item{mCn}{
  	mCn, object created in patchwork.CG.copynumbers().
	}
	\item{t}{
  	A list between int and ai. (ratio and allelic imbalance)
	}
	\item{name}{
  	Default is ''. When method is run by patchwork.CG.copynumbers() it will inherit that functions name parameter.
	}
  	\item{xlim}{
  	Default is c(0,2.4). Limit of x axis.
	}
  	\item{ylim}{
  	Default is c(0,1). Limit of y axis.
	}
}


\details{
Vertical axis: Allelic Imbalance. \cr
Horizontal axis: Relative coverage. \cr

The naming scheme is Copynumber-m-LesserAlleleDistribution
So for example 2m0 means copynumber = 2, both allels are the same
whereas 2m1 means copynumber = 2, 1 allele each.

Another example:
4m0, copynumber = 4, All allels are the same. (Loss of heterozygosity).
4m1, copynumber = 4, 3 alleles are the same, one is different.
4m2, copynumber = 4, 2 alleles each.

The total number of alleles present are always the copynumber.
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
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.CG.copynumbers}}
}

%%\examples{

%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line