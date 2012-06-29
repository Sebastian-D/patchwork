 \name{patchwork.help}
\alias{patchwork.readme}
\alias{readme.patchwork}
\alias{patchwork.help}
\alias{help.patchwork}


\title{
Welcome to patchwork
}

\description{

patchwork allows you to obtain allele-specific copy number information from BAM files.

To use patchwork on CompleteGenomics data download patchworkCG.

//////////////////////////////////////////////// \cr
http://patchwork.r-forge.r-project.org/ \cr
////////////////////////////////////////////////

It is highly recommended that you visit the webpage for this project as it contains the most recently updated
information and instructions regarding everything patchwork.

}

%\usage{
%patchwork.plot(BamFile,Pileup,reference=NULL,normal.bam=NULL,normal.pileup=NULL,Alpha=0.0001,SD=1)
%}

%\arguments{
%  \item{BamFile}{
%  	Aligned, Sorted Bam-file.
% 	}
%   \item{Pileup}{
% 	Pileup file generated through 
% 	samtools -vcf reference.fasta bamfile > outfile.
% 	}
%   \item{reference}{
%   	Default is NULL.
% 	Path to a reference file that can be created using patchwork.createreference() or downloaded from patchworks homepage.
% 	}
% 	\item{normal.bam}{
% 	Default is NULL.
% 	The matched normal sample of the your BamFile.
% 	}
% 	\item{normal.pileup}{
% 	Default is NULL.
% 	The pileup of your normal sample.
% 	}
% 	\item{Alpha}{
% 	Default 0.0001, change if you want to try to get a better segmentation from patchwork.segment().
% 	From DNAcopy (?segment):
% 	 alpha: significance levels for the test to accept change-points.
% 	}
% 	\item{SD}{
% 	Default 1, change if you want to try to get a better segmentation from patchwork.segment().
% 	From DNAcopy (?segment):
% 	undo.SD: the number of SDs between means to keep a split if undo.splits="sdundo".
% 	}
% }

\details{

//////////////////////////// \cr
Installation Guide: \cr
////////////////////////////

For patchwork to run correctly, or even install, you will need to install DNACopy from Bioconductor.

Start R (as of writing version 2.14.2)

Text with a ">" infront is R executions.

    > source("http://bioconductor.org/biocLite.R") \cr
    > biocLite("DNAcopy")

After DNAcopy is installed you can install patchwork and patchworkData using these commands:

    > install.packages("patchworkData", repos="http://R-Forge.R-project.org") \cr
    > install.packages("patchwork", repos="http://R-Forge.R-project.org")

If for some reason that does not work add the 'type="source"' to it, as so:

    > install.packages("patchworkData", repos="http://R-Forge.R-project.org",type="source") \cr
    > install.packages("patchwork", repos="http://R-Forge.R-project.org",type="source")

The source codes can also be downloaded outside of R at https://r-forge.r-project.org/R/?group_id=1250 

//////////////////////////// \cr
	Tutorial: \cr
////////////////////////////

	///////////////////// \cr
		Requirements: \cr
	/////////////////////

-A working R installation.

-Your samples sorted BAM and BAI files: \cr
	To get the BAI file run "samtools index <yourfile>.bam" \cr
	To sort a bamfile "samtools sort <yourfile>.bam <sortedfile>.bam"

-A pileup of your sample. \cr
To do this step you will need a reference human genome fasta file, for example HumanGenome19.fasta.

At the time of writing there should be several places you can find this; The Sanger institute, 1000genomes project etc. You probably already have one as you needed one previously for the alignment of your tumor sequence.

Use SAMtools, version 0.1.18 or older, to produce the pileup: \cr
	"samtools pileup -vcf <reference.fasta> <yourfile>.bam > <yourpileup>"

The pileup should have this format:

less pileup \cr \tabular{rlllllllll}{
	chr1 \tab 10179 \tab   c    \tab    W    \tab    0    \tab    0   \tab     60   \tab    5    \tab    ,,tna  \tab  ##6!3 \cr
	chr1  \tab   10180  \tab  t  \tab      W     \tab   6   \tab     6   \tab     60  \tab     5   \tab     ,,,aa \tab   ##369 \cr
	chr1  \tab   10377  \tab  a   \tab     R    \tab    0   \tab     1   \tab     60  \tab     1   \tab     g   \tab     @ \cr
	chr1  \tab   11391  \tab  t   \tab     A    \tab    0   \tab     3   \tab     60  \tab     1   \tab     a   \tab     B \cr
	chr1 \tab    18592  \tab  C   \tab     Y     \tab   0   \tab     3   \tab     60   \tab    1   \tab     t   \tab     B \cr
	chr1   \tab  23359  \tab  c   \tab     S \tab       0   \tab     2   \tab     60  \tab     1   \tab     g   \tab     A \cr
	chr1 \tab    24067  \tab  A   \tab     R      \tab  0   \tab     2    \tab    60   \tab    1   \tab     G   \tab     A \cr
	chr1  \tab   30315  \tab  G  \tab      C   \tab     3   \tab     24   \tab    60  \tab     2  \tab      cC  \tab     =; \cr
	chr1  \tab   92200  \tab  a   \tab     W   \tab    0    \tab    3    \tab    60  \tab     1    \tab    t   \tab     B \cr
	chr1  \tab   96592  \tab  T   \tab     C   \tab     0   \tab     3    \tab    60  \tab     1   \tab     C  \tab     B \cr
	chr1  \tab   100140 \tab  a   \tab     M    \tab    0   \tab     1    \tab    60  \tab     1   \tab     C  \tab      @ \cr
	chr1  \tab   104697 \tab  g   \tab     K    \tab    0   \tab     2    \tab    60  \tab     1    \tab    T  \tab      A \cr
	chr1  \tab   127285 \tab  A  \tab      R   \tab     0   \tab     3   \tab     60  \tab     1   \tab     g  \tab      B \cr
}

-(optional) A matched normal sample to your tumor in BAM \cr
-(optional) A pileup of your normal sample BAM \cr
-(optional) A standard Reference file. (Illumina/Solexa, SOLiD or your own)

What we mean by this (optional) tag is that you do not need to have all these
files. There is really no point in using all of them. 

If you have a matched normal sample you should make a pileup of this and use those two arguments.

In some cases it may be better to use the pileup and a reference file.

You can also run patchwork.plot with only a reference file.
See patchwork.createreference for information about the reference file.

	/////////////////////////////////////////////// \cr
		Execution: patchwork.createreference() \cr
	///////////////////////////////////////////////

This function creates a reference using a pool of samples of your selection. Choosing the samples to use for reference creation you should consider these factors: \cr
	- They should be sequenced using the same technique \cr
	- They should be from the same organism \cr
	- They should be non-tumerous

It is recommended that you use at least 3 bam files to create your own reference.

Start R

Load the patchwork and patchworkData libraries

	> library(patchwork) \cr
	> library(patchworkData)

Read the amazing documentation for patchwork.createreference()
	
	> ?patchwork.createreference

Execute the function, pointing to your desired files.

	>patchwork.createreference("file1","path/to/file2","../file3","~/heres/file4",output="REFOUT")

This will generate REFOUT.Rdata, or whichever prefix you chose, in your working directory. Use this file for the reference argument of patchwork.plot().

	//////////////////////////////////// \cr
		Python and Pysam \cr
	////////////////////////////////////

Patchwork uses a python mod called pysam for chromosome reading.
The package already checks if you have it installed and otherwise installs
it for you the first time you run patchwork however if you lack python
or python development on your system the installation fails. It seems to already
be included in MAC OS X but on other systems you should install/update it.

See http://wiki.python.org/moin/BeginnersGuide/Download

	//////////////////////////////////////// \cr
		Execution: patchwork.plot() \cr
	////////////////////////////////////////

It is recommended that you run patchwork from a "clean" working directory. In this way you do not run the risk of having files write over eachother. If you do not want to type paths into R you may also want to put the required files in this folder.

Execution may take quite a while depending on the size of your sample! So if possible run it on a dedicated computer.


Start R

Load the patchwork and patchworkData libraries

	> library(patchwork) \cr
	> library(patchworkData)

Read the excellent documentation for patchwork.plot()

	> ?patchwork.plot

Excerpt from ?patchwork.plot:

////////////////////////////////////////\cr
	Usage:

	     patchwork.plot(BamFile,Pileup,reference=NULL,normal.bam=NULL,normal.pileup=NULL,Alpha=0.0001,SD=1)
	     
	Arguments:

		 BamFile: Aligned, Sorted Bam-file.

		  Pileup: Pileup file generated through samtools -vcf reference.fasta
		          bamfile > outfile.

		reference: Default is NULL.  Path to a reference file that can be
		          created using patchwork.createreference().

		normal.bam: Default is NULL.  The matched normal sample of the your
		          BamFile.

		normal.pileup: Default is NULL.  The pileup of your normal sample.

		   Alpha: Default 0.0001, change if you want to try to get a better
		          segmentation from patchwork.segment().  From DNAcopy
		          (?segment): alpha: significance levels for the test to accept
		          change-points.

		      SD: Default 1, change if you want to try to get a better
		          segmentation from patchwork.segment().  From DNAcopy
		          (?segment): undo.SD: the number of SDs between means to keep
		          a split if undo.splits="sdundo".
//////////////////////////////////////// \cr


Perform patchwork.plot() with desired arguments.

	>patchwork.plot(BamFile="patchwork.example.bam",Pileup="patchwork.example.pileup",reference="../HCC1954/datasolexa.RData") \cr
	Initiating Allele Data Generation\cr
	Initiating Read Chromosomal Coverage \cr
	Reading chr1 \cr
	Reading chr2 \cr
	Reading chr3 \cr
	Reading chr4 \cr
	Reading chr5 \cr
	Reading chr6 \cr
	Reading chr7 \cr
	Reading chr8 \cr
	Reading chr9 \cr
	Reading chrX \cr
	Reading chrY \cr
	Reading chr10 \cr
	Reading chr11 \cr
	Reading chr12 \cr
	Reading chr13 \cr
	Reading chr14 \cr
	Reading chr15 \cr
	Reading chr16 \cr
	Reading chr17 \cr
	Reading chr18 \cr
	Reading chr19 \cr
	Reading chr20 \cr
	Reading chr21 \cr
	Reading chr22 \cr
	Read Chromosomal Coverage Complete \cr
	Initiating GC Content Normalization \cr
	GC Content Normalization Complete \cr
	Initiating Smoothing \cr
	Smoothing Chromosome: chr1 \cr
	Smoothing Chromosome: chr2 \cr
	Smoothing Chromosome: chr3 \cr
	Smoothing Chromosome: chr4 \cr
	Smoothing Chromosome: chr5 \cr
	Smoothing Chromosome: chr6 \cr
	Smoothing Chromosome: chr7 \cr
	Smoothing Chromosome: chr8 \cr
	Smoothing Chromosome: chr9 \cr
	Smoothing Chromosome: chrX \cr
	Smoothing Chromosome: chrY \cr
	Smoothing Chromosome: chr10 \cr
	Smoothing Chromosome: chr11 \cr
	Smoothing Chromosome: chr12 \cr
	Smoothing Chromosome: chr13 \cr
	Smoothing Chromosome: chr14 \cr
	Smoothing Chromosome: chr15 \cr
	Smoothing Chromosome: chr16 \cr
	Smoothing Chromosome: chr17 \cr
	Smoothing Chromosome: chr18 \cr
	Smoothing Chromosome: chr19 \cr
	Smoothing Chromosome: chr20 \cr
	Smoothing Chromosome: chr21 \cr
	Smoothing Chromosome: chr22 \cr
	Smoothing Complete \cr
	Initiating Segmentation \cr 
	Note: If segmentation fails to initiate the probable reason is that you have not installed the R package DNAcopy. See the homepage, http://patchwork.r-forge.r-project.org/ , or ?patchwork.readme for installation instructions.  \cr
	Analyzing: chr1.p \cr
	Analyzing: chr1.q \cr
	Analyzing: chr2.p \cr
	Analyzing: chr2.q \cr
	Analyzing: chr3.p \cr
	Analyzing: chr3.q \cr
	Analyzing: chr4.p \cr
	Analyzing: chr4.q \cr
	Analyzing: chr5.p \cr
	Analyzing: chr5.q \cr
	Analyzing: chr6.p \cr
	Analyzing: chr6.q \cr
	Analyzing: chr7.p \cr
	Analyzing: chr7.q \cr
	Analyzing: chr8.p \cr
	Analyzing: chr8.q \cr
	Analyzing: chr9.p \cr
	Analyzing: chr9.q \cr
	Analyzing: chrX.p \cr
	Analyzing: chrX.q \cr
	Analyzing: chrY.p \cr
	Analyzing: chrY.q \cr
	Analyzing: chr10.p \cr
	Analyzing: chr10.q \cr
	Analyzing: chr11.p \cr
	Analyzing: chr11.q \cr
	Analyzing: chr12.p \cr
	Analyzing: chr12.q \cr
	Analyzing: chr13.q \cr
	Analyzing: chr14.q \cr
	Analyzing: chr15.q \cr
	Analyzing: chr16.p \cr
	Analyzing: chr16.q \cr
	Analyzing: chr17.p \cr
	Analyzing: chr17.q \cr
	Analyzing: chr18.p \cr
	Analyzing: chr18.q \cr
	Analyzing: chr19.p \cr
	Analyzing: chr19.q \cr
	Analyzing: chr20.p \cr
	Analyzing: chr20.q \cr
	Analyzing: chr21.p \cr
	Analyzing: chr21.q \cr
	Analyzing: chr22.q \cr
	Segmentation Complete \cr
	Initiating Segment data extraction (Medians and AI) \cr
	Segment data extraction Complete \cr
	 
	 
	Saving information objects needed for patchwork.copynumbers in copynumbers.Rdata 
	 
	 
	Initiating Plotting \cr
	Plotting Complete \cr
	Shutting down..... \cr
	Warning messages:


If you did it correctly it should have generated similar output as can be seen above.
Your working directory should now have the plots generated from the function, 1 overhead plot and 24 chromosomal plots.
The working directory should also contain the files:\cr
	- copynumbers.Rdata\cr
	- data.Rdata\cr
	- pile.alleles\cr
	- pile.alleles.Rdata\cr
	- Segments.Rdata\cr
	- smoothed.Rdata\cr

These were created for swifter re-runs of the function should something unforseen happen during execution. For example if something goes wrong during a final step it would be a terrible hassle to read all the chromosomes again.

!!!THIS DOES MEAN THAT TO COMPLETELY RE-RUN PATCHWORK FROM SCRATCH IN THE SAME WORKING DIRECTORY ALL OF THE ABOVE FILES NEED TO BE DELETED!!!

For some details of the plot functions see their documentation:

	> ?karyotype \cr
	> ?karyotype_chroms

To interpret the plots here is an excerpt from their documentation:

	///////////////////////////\cr
		karyotype()\cr
	///////////////////////////

	Description:

	     Plots each,color coded by chromosomal coordinate, chromosome
	     against a background of the complete genome.

    Details:

	     Vertical axis: Allelic imbalance. \cr
	     Horizontal axis: Total intensity. \cr
	     The plot is a overview, for a closer look see the plots generated
	     by karyotype_chroms().
	
	///////////////////////////\cr
		karyotype_chroms()\cr
	///////////////////////////

	Description:

	     Visualises the calculated data of patchwork.plot() for each
	     chromosome. 

    Details:

	     *TOP* \cr
	     Vertical axis: Allelic Imbalance
	     Horizontal axis: Total Intensity
	     The chromosome plotted against the complete genome background. The
	     separation between clusters within the plot are due to the
	     fluctuating intensity and allelic imbalance and as such display
	     the varying allele counts and copy numbers. Longer/larger segments
	     have bigger circles. Darker circles show more content as they are
	     ontop of eachother.

	     *MIDDLE* \cr
	     Vertical axis: Total Intensity
	     Horizontal axis: Chromosomal coordinate
	     The chromosome in questions total intensity plotted against the
	     position on the chromosome.

	     *BOTTOM* \cr
	     Vertical axis: Allelic Imbalance
	     Horizontal axis: Chromosomal coordinate
	     The chromosome in questions allelic imbalance plotted against the
	     position on the chromosome.


	//////////////////////////////////////////////\cr
		Execution: patchwork.copynumbers()\cr
	//////////////////////////////////////////////

The only file you absolutely must have in your working for the next part of execution is copynumbers.Rdata.

Read the outstanding documentation for patchwork.copynumbers()

	> ?patchwork.copynumbers

The cliffnotes are that you will need to look at the plots generated by patchwork.plot() to be able to correctly input the
arguments for patchwork.copynumbers().

Excerpt from patchwork.copynumbers() documentation:

//////////////////////////////////////////////\cr
	Usage:

	     patchwork.copynumbers(cn2,delta,het,hom,maxCn=8,ceiling=1,forcedelta=F)
	     
	Arguments:

	    name: Default is "copynumbers_". First part of output name for plots
	          generated from patchwork.copynumbers().

	    cn2: The approximate position of copy number 2,diploid, on total
	          intensity axis.

	   	delta: The difference in total intensity between consecutive copy
	          numbers. For example 1 and 2 or 2 and 3.  If copy number 2
	          has total intensity ~0.6 and copy number 3 har total
	          intensity ~0.8 then delta would be 0.2.

	    het: Allelic imbalance ratio of heterozygous copy number 2.

	    hom: Allelic imbalance ratio of Loss-of-heterozygosity copy number
	          2.

	   	maxCn: Highest copy number to calculate for. Default is 8.

	 	ceiling: Default is 1.

	 	forcedelta: Default is FALSE. If TRUE the delta value will not be
          subject to small adjustment changes.
//////////////////////////////////////////////


On the projects webpage, http://patchwork.r-forge.r-project.org/  (patchwork tab -> execution tab), there will also be a picture accompanying this portion pointing to the argument values should the documentation for the function not be sufficient.

	>patchwork.copynumbers(name="Example_",cn2=0.8,delta=0.3,het=0.5,hom=0.8)

patchwork.copynumbers() generates plots and should not take a large amount of time to complete. When it is finished you should have 24 chromosome plots, one for each chromosome, and an overview plot with the approximate positions of copy numbers and allele ratios.

For some details of the plot functions see their documentation:

	> ?karyotype_check \cr
	> ?karyotype_chromsCN

To interpret the plots here is an excerpt from their documentation:
	
///////////////////////////////////\cr
		karyotype_check()\cr
///////////////////////////////////

	Description:

	     Plots the whole genome coverage vs allelic imbalance with the
	     approximated areas copynumbers and allele constitution tagged.

	Details:

	     Vertical axis: Allelic Imbalance \cr
	     Horizontal axis: Relative coverage

	     The naming scheme is Copynumber-m-LesserAlleleDistribution. \cr
	     For example 2m0 means copynumber = 2, both alleles are the same whereas
	     2m1 means copynumber = 2, 1 allele each.

	     Another example: 4m0, copynumber = 4, All allels are the same.
	     (Loss of heterozygosity). 4m1, copynumber = 4, 3 alleles are the
	     same, one is different. 4m2, copynumber = 4, 2 alleles each.

	     The total number of alleles present are always the copynumber.

/////////////////////////////////// \cr
		karyotype_chromsCN() \cr
///////////////////////////////////

	Description:

	     Visualises the calculated data of patchwork.plot() +
	     patchwork.copynumbers() for each chromosome. See details for a
	     walkthrough of the plot.

	Details:

	     *TOP* \cr
	     Vertical axis: Allelic Imbalance
	     Horizontal axis: Total Intensity
	     The chromosome plotted against the complete genome background. The
	     separation between clusters within the plot are due to the
	     fluctuating intensity and allelic imbalance and as such display
	     the varying allele counts and copy numbers. Longer/larger segments
	     have bigger circles. Darker circles show more content as they are
	     ontop of eachother.

	     *TOP MIDDLE* \cr
	     Vertical axis: Copynumber
	     Horizontal axis: Chromosomal coordinate
	     Displays the total and minor copynumbers for different segments of
	     the chromosome in question.

	     *LOWER MIDDLE* \cr
	     Vertical axis: Total Intensity
	     Horizontal axis: Chromosomal coordinate
	     The chromosome in questions total intensity plotted against the
	     position on the chromosome.

	     *BOTTOM* \cr
	     Vertical axis: Allelic Imbalance
	     Horizontal axis: Chromosomal coordinate
	     The chromosome in questions allelic imbalance plotted against the
	     position on the chromosome.



If you have any questions please feel free to contact us and we will help you to the best of our extent!

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
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se} \cr
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

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