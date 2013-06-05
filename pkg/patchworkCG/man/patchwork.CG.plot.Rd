\name{patchwork.CG.plot}
\alias{patchworkCG}
\alias{patchwork.CG.plot}
\alias{patchwork.CG}


\title{
The main function of patchworkCG.
}

\description{
Reads files and performs calculations later used in the plotting functions.

In details there is some very important information, highly recommended!
}

\usage{
patchwork.CG.plot(path=NULL,name='CG_sample',manual_file_input=FALSE,
			masterVarBeta=NULL,somaticCnvSegments=NULL,
			depthOfCoverage=NULL)
}

\arguments{
  \item{path}{
  	Path to the ASM folder. One of the subfolders of your completegenomics directory. Should be given as a string.
    Ex. path="path/to/ASM".
	}
  \item{name}{
	Default is 'CG_sample'. The name you wish associated with the plots that will be generated.
	}
  \item{manual_file_input}{
	Default is FALSE. If you set it to TRUE you will be prompted to provide path and filename for
	the masterVarBeta, somaticCnvSegmentsNondiploid and depthOfCoverage file.
	}
  \item{masterVarBeta}{
  	Path to and COMPLETE NAME of the masterVarBeta file, should you wish to implement it directly. Useful if you want to run
  	the program on multiple samples and wish to create a script going from file to file. This saves you
  	from having to conserve the CompleteGenomics file structure.
  	Default is NULL.
  	}
  \item{somaticCnvSegments}{
  	Path to and COMPLETE NAME of the somaticCnvSegments file, should you wish to implement it directly. Useful if you want to run
  	the program on multiple samples and wish to create a script going from file to file. This saves you
  	from having to conserve the CompleteGenomics file structure.
  	Default is NULL.
  	}
  	\item{depthOfCoverage}{
  	Path to and COMPLETE NAME of the depthOfCoverage file, should you wish to implement it directly. Useful if you want to run
  	the program on multiple samples and wish to create a script going from file to file. This saves you
  	from having to conserve the CompleteGenomics file structure.
  	Default is NULL.
  	}
}

\details{
Before execution make sure that you have unpacked the necessary files. masterVarBeta, located in your ASM/
folder,somaticCnvSegmentsNonDiploidBeta, located in your ASM/CNV/ folder, and depthOfCoverage, also located in your ASM/CNV/ folder.

If you wish to use this function in a script so it performs its function on multiple genome entries, remember
to also loop over the "name" parameter or you run the risk of the plots replacing eachother. The easiest way
to script this to run over multiple files is to use the masterVarBeta, somaticCnvSegments and depthOfCoverage
parameters as opposed to the path parameter.

It is prudent not do run multiple instances of patchwork.CG.plot in the same directory. This can cause issues
with the generated file, read below, and also with plot names if you have kept the default value.

During execution of patchwork.CG.plot a file named CG.Rdata will be built then saved in your working directory. This file is used
for plotting purposes in patchwork.CG.copynumbers(). As such, patchwork.CG.copynumbers() should either be executed in the same working
directory as this file is located or using the CGfile parameter of patchwork.CG.copynumbers().
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{

#Load patchworkCG
>library(patchworkCG)

#An example of running patchwork.CG.plot. Note that the path is to the ASM folder
#of the default CompleteGenomics File structure.
>patchwork.CG.plot(path="Data/HCC2218/GS00258-DNA_B03/GS00258-DNA_B03/ASM",
				name="HCC2218_")

Reading files from ASM folder 
File input complete 
Performing calculations 
Calculations complete 
Saving objects to CG.Rdata 
Initiating Plotting 
Plotting Complete 
patchwork.CG.plot Complete 
Please read documentation on running patchwork.CG.copynumbers

#The values used to generate the plots are saved in three objects,
#mastervar,segs and depcov, and can be viewed loading CG.Rdata.
>load("CG.Rdata")

> head(mastervar)
      chr begin   end vartype ref_count tot_count ref_countN tot_countN
876  chr1 46669 46670     snp         1         1          2          3
912  chr1 49313 49315     sub        69       129         71        123
1030 chr1 55387 55388     snp        64       211         78        247
1042 chr1 55815 55816     snp         5        12          3         10
1054 chr1 56484 56485     snp         0        11          2          6
1082 chr1 57989 57990     snp        58       187         56        177
     mut_count max min      ratio
876          0   1   0 0.01635786
912         60  69  60 2.11016381
1030       147 147  64 3.45150825
1042         7   7   5 0.19629431
1054        11  11   0 0.17993645
1082       129 129  58 3.05891963

> head(segs)
   chr   start     end avgnormcov relcov     ratio snvs        ai
1 chr1   10000  177417       88.2   1.36 1.0434479   42 0.5449591
2 chr1  227417  267719       88.2   1.36 1.0434479    0 0.0000000
3 chr1  317719  471368       88.2   1.36 1.0434479    2 1.0000000
4 chr1  521368  600000       61.4   0.95 0.8015226    1 1.0000000
5 chr1  600000 2634220       93.9   1.45 1.0965535  895 0.6904948
6 chr1 2684220 3845268       93.9   1.45 1.0965535  570 0.8107859

> head(depcov)
   chr  begin    end avgnormcov     ratio
1 chr1  10000 100000       94.8 1.0897044
2 chr1 100000 177417      102.9 1.1622234
3 chr1 227417 267719      124.4 1.3547120
4 chr1 317719 400000       99.5 1.1317833
5 chr1 400000 471368       83.9 0.9921172
6 chr1 521368 600000       70.1 0.8685663 
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line