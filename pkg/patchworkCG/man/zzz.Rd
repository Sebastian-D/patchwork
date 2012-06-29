\name{patchworkCG.readme}
\alias{patchworkCG.readme}
\alias{readme.patchworkCG}
\alias{patchworkCG.help}
\alias{help.patchworkCG}

\title{
Welcome to patchworkCG
}

\description{
PatchworkCG allows you to obtain allele-specific copy number information from Complete Genomics data.

To use patchwork on other data formats,bam and pileups, download patchwork and patchworkData.

//////////////////////////////////////////////// \cr
http://patchwork.r-forge.r-project.org/ \cr
////////////////////////////////////////////////

It is highly recommended that you visit the webpage for this project as it contains the most recently updated information and instructions regarding everything patchwork.
}



\details{

//////////////////////////////////////////////// \cr
Installation Guide: PatchworkCG \cr
//////////////////////////////////////////////// 

Start the latest version of R.(as of writing 2.14.2)

Text with a ">" infront is R executions.

	> install.packages("patchworkCG", repos="http://R-Forge.R-project.org")

If for some reason that does not work add the 'type="source"' to it, as so:

	> install.packages("patchworkCG", repos="http://R-Forge.R-project.org",type="source")

The source codes can also be downloaded outside of R at https://r-forge.r-project.org/R/?group_id=1250 

////////////////// \cr
	Tutorial: \cr
//////////////////

	////////////////// \cr
		Requirements: \cr
	////////////////// 

An unpacked CompleteGenomics tumor genome data ASM directory. Specifically the files ASM/masterVarBeta, CNV/depthOfCoverage and CNV/somaticCnvSegmentsNondiploidBeta.

Some helpful commands for unpacking your CompleteGenomics data (UNIX): \cr
To unpack .tar files:tar -zxvf filename.tar \cr
To unpack .bz2 files:bunzip2 filename.bz2

////////////////////////////////// \cr
Execution: patchwork.CG.plot() \cr
//////////////////////////////////


Start R in your working directory. (Or use setwd() after starting R wherever you choose)

	> library(patchworkCG)

See 

	> ?patchwork.CG.plot

Excerpt from ?patchwork.CG.plot:

//////////////////////////////////////////// \cr	
	Usage:

     	patchwork.CG.plot(path=NULL,name='CG_sample',manual_file_input=FALSE,
                             masterVarBeta=NULL,somaticCnvSegments=NULL,
                             depthOfCoverage=NULL)
     
	Arguments:

	    path: Default is NULL. Path to the ASM folder. One of the subfolders of
	    	your completegenomics directory.

	    name: Default is 'CG_sample'. The name you wish associated with
	          the plots that will be generated.

		manual_file_input: Default is FALSE. If you set it to TRUE you will be
	          prompted to provide path and filename for the masterVarBeta,
	          somaticCnvSegmentsNondiploid and depthOfCoverage file.

		masterVarBeta: Path to AND COMPLETE name of the masterVarBeta file,
	          should you wish to implement it directly. Useful if you want
	          to run the program on multiple samples and wish to create a
	          script going from file to file. This saves you from having to
	          conserve the CompleteGenomics file structure.  Default is
	          NULL.

		somaticCnvSegments: Path to AND COMPLETE name of the somaticCnvSegments
	          file, should you wish to implement it directly. Useful if you
	          want to run the program on multiple samples and wish to
	          create a script going from file to file. This saves you from
	          having to conserve the CompleteGenomics file structure.
	          Default is NULL.

		depthOfCoverage: Path to AND COMPLETE name of the depthOfCoverage file,
	          should you wish to implement it directly. Useful if you want
	          to run the program on multiple samples and wish to create a
	          script going from file to file. This saves you from having to
	          conserve the CompleteGenomics file structure.  Default is
	          NULL.
////////////////////////////////////////////

Probably the easiest way of running patchwork.CG.plot is pointing to the ASM folder:

	> patchwork.CG.plot(path="path/to/ASM/") \cr
	Reading files from ASM folder \cr
	File input complete \cr
	Performing calculations \cr
	Calculations complete \cr
	Saving objects to CG.Rdata \cr
	Initiating Plotting \cr
	Plotting Complete \cr
	patchwork.CG.plot Complete \cr
	 Please read documentation on running patchwork.CG.copynumbers 

Once it is complete there should be several plots in your working directory as
well as the file "CG.Rdata".

Execution should take no longer half an hour.

	//////////////////////////////////////////// \cr
	Execution: patchwork.CG.copynumbers() \cr
	////////////////////////////////////////////

For this function to work you will need to have obtained the CG.Rdata file from
patchwork.CG.plot(). If it is located in your working directory you will not need to
point to it, otherwise you can input it with the argument "CGfile".

See 

	> ?patchwork.CG.copynumbers 

Excerpt from ?patchwork.CG.copynumbers:

//////////////////////////////////////////// \cr
	Usage:

	    patchwork.CG.copynumbers(cn2,delta,het,hom,maxCn=8,ceiling=1,name="copynumbers_",
	                                     CGfile=NULL,forcedelta=F)
     
	Arguments:

	    cn2: The approximate position of copy number 2, diploid, on total
	          intensity axis.

	   	delta: The difference in total intensity between consecutive copy
	          numbers. For example 1 and 2 or 2 and 3.  If copy number 2
	          has total intensity ~0.6 and copy number 3 har total
	          intensity ~0.8 then delta would be 0.2.

	    het: Allelic imbalance ratio of heterozygous copy number 2.

	    hom: Allelic imbalance ratio of Loss-of-Heterozygosity copy number
	          2.

	   	maxCn: Highest copy number to calculate for. Default is 8.

	 	ceiling: Default is 1.

	    name: Default is "copynumbers_". The name you want attached to generated
	          plots.

	  	CGfile: Default is NULL. If your CG.Rdata file is not in your working
	          directory, and you dont wish to move it to your working
	          directory, you can simply input the path here as CGfile =
	          "path/to/file/CG.Rdata" so patchwork.CG.copynumbers() can find its
	          data.

	    forcedelta: Default is FALSE. If TRUE the delta value will be absolute
          	and not subject to adjustments.
////////////////////////////////////////////

On our webpage there will also be a picture accompanying this portion pointing to the argument values should the documentation for the function not be sufficient.
(http://patchwork.r-forge.r-project.org/)

	> patchwork.CG.copynumbers(cn2=0.8, delta=0.3, het=0.2, hom=0.8, name="Example_", CGfile="path/to/CG.Rdata")

This should have generated several plots in your working directory.


	///////////////////////////////// \cr
		Interpreting plots \cr
	/////////////////////////////////

It is recommended that you visit our project homepage for guidance in
interpreting the plots you have generated using patchworkCG.

	http://patchwork.r-forge.r-project.org/

Other aids in this is the documentation
	
	>?CG_Ka_check \cr
	>?CG_KaCh \cr
	>?CG_KaChCN \cr
	>?CG_karyotype


If you have any questions please feel free to contact us and we will help you to the best of our extent!

}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se} \cr
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\examples{
%\dontrun{}
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%Run examples of karyotype stuff. or not, maybe just point back to CNA.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line