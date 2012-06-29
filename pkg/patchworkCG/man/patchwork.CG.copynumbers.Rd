\name{patchwork.CG.copynumbers}
\alias{patchwork.CG.copynumbers}
\alias{copynumbers}


\title{
Function plot allele-specific copy numbers of the genome
}

%Markus:briefly what is done in this function?
\description{
Input data for this function should be assessed from any of the CG_KaCh()
plots.

See argument details.

\bold{Very important: For this function to work you must have the file CG.Rdata,
generated using patchwork.CG.plot(), in your working directory!}
}

\usage{
patchwork.CG.copynumbers(cn2,delta,het,hom,maxCn=8,ceiling=1,name="copynumbers_",
				CGfile=NULL,forcedelta=F)
}

\arguments{
  \item{cn2}{
	The approximate position of copy number 2, diploid, on total intensity axis.
	}
  \item{delta}{
	The difference in total intensity between consecutive copy numbers. For example 1 and 2 or 2 and 3.
	If copy number 2 has total intensity ~0.6 and copy number 3 har total intensity ~0.8 then delta would be 0.2.
	}
  \item{het}{
	Allelic imbalance ratio of heterozygous copy number 2.
	}
%What if LOH cn2 and het cn2 isnt in the sample?
  \item{hom}{
	Allelic imbalance ratio of Loss-of-Heterozygosity copy number 2.
	}
  \item{maxCn}{
	Highest copy number to calculate for. Default is 8.
	}
  \item{ceiling}{
	Default is 1.
	}
  \item{name}{
	Default is "copynumbers_". The name you want attached to generated plots.
	}
	\item{CGfile}{
	Default is NULL. If your CG.Rdata file is not in your working directory, and you dont wish to move
	it to your working directory, you can simply input the path here as CGfile = "path/to/file/CG.Rdata"
	so patchwork.CG.copynumbers() can find its data.
	}
	\item{forcedelta}{
	Default is FALSE. If TRUE the delta value will be absolute and not subject to adjustments.
	}
}

\details{
1. Run the patchwork.CG.plot correctly. This will generate CG.Rdata. \cr
2. View your plots and decide the correct parameters for the other arguments of this function. \cr
3. Make sure CG.Rdata is in your working directory, or set the CGfile parameter to point to CG.Rdata. \cr
4. Run the function. \cr

Upon completion two files will be generated, copynumbers.Rdata which contains the objects used to generate the plots and Copynumbers.csv.
Other than this plots from CG_Ka_check and CG_KaChCN will be generated.
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
\code{\link{patchwork.CG.plot}} \cr
\code{\link{subfunctions_copynumbers}} \cr
\code{\link{CG_Ka_check}} \cr
\code{\link{CG_KaChCN}}
}

\examples{
\dontrun{
>library(patchworkCG)

#With CG.Rdata in working directory. Values for parameters estimated from
#HCC2218 chromosome 1. 
#See plot patchworkCG/example_plots/HCC2218_KaCh_chr1.png
>patchwork.CG.copynumbers(cn2=0.85,delta=0.35,het=0.25,hom=0.95,name="HCC2218_")

#Pointing to location of CG.Rdata file
>patchwork.CG.copynumbers(cn2=0.85,delta=0.35,het=0.25,hom=0.95,name="HCC2218_",
				CGfile="path/to/CG.Rdata")

#Load the objects generated from patchwork.CG.copynumbers
#For comparison purposes only
>load("copynumbers.Rdata")

> head(ai)
$cn2m1
[1] 0.2460513

$cn2m0
[1] 0.95

$cn1m0
[1] 0.9690196

$cn0m0
[1] NA

$cn3m1
[1] 0.2677295

$cn3m0
[1] 0.9671712

> head(int)
$cn2
[1] 0.85

$cn1
[1] 0.5

$cn0
[1] 0.15

$cn3
[1] 1.2

$cn4
[1] 1.55

$cn5
[1] 1.9

> head(regions)
   chr  start    end avgnormcov relcov     ratio snvs        ai np Cn mCn
1 chr1  10000 177417       79.0   1.05 0.9660335   24 0.5531915 17  2   1
2 chr1 227417 267719       79.0   1.05 0.9660335    2 0.3243243  4  2   1
3 chr1 317719 471368       79.0   1.05 0.9660335    2 0.7826087 15  2   0
4 chr1 521368 600000      170.8   2.26 2.0792721    2 0.1944444  8  6   3
5 chr1 600000 700000       74.6   0.99 0.9108316    2 0.5333333 10  2   1
6 chr1 700000 900000       51.9   0.69 0.6348220  109 0.4960325 20  1   0
  fullCN
1  cn2m1
2  cn2m1
3  cn2m0
4  cn6m3
5  cn2m1
6  cn1m0
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line