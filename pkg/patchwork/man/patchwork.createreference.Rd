\name{patchwork.createreference}
\alias{createreference}
\alias{patchwork.createreference}


\title{
Patchwork function to create your own normal reference file.
}

\description{
Create your own normal reference file for use with patchwork.plot(). This may have advantages over using the
standard issue references. It is recommended to use atleast 3 bamfiles to create a
reference file.
}

\usage{
patchwork.createreference(\dots,output="REFOUT")
}

\arguments{
  \item{\dots}{
  	The full path and filename to all of the BAM files
  	that you wish to use for creating the reference file.
  	ex patchwork.createreference("file1.bam","file2.bam","file3.bam",etc)
	}
	\item{output}{
	Default is "REFOUT". This will result in REFOUT.RData being created in your working
	directory.
	}
}

\details{
This function will take as its arguments a number of bamfiles that you
specify, normalize them according to patchwork.readChroms and
patchwork.GCNorm. REFOUT.Rdata, or whichever name you chose, is created in your
working directory for further use. This file should then
be sent in patchwork.plot to be used as a reference file.

patchwork.plot(bamfile,pileup,reference="path/to/your/REFOUT.Rdata")
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

%\code{\link{help}}
\seealso{
\code{\link{patchwork.plot}} \cr
\code{\link{patchwork.readChroms}} \cr
\code{\link{patchwork.applyref}}	\cr
\code{\link{patchwork.GCNorm}}
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