\name{TAPS_plot}
\alias{TAPS_plot}
\alias{Taps_plot}
\alias{TAPS.plot}
\alias{Taps.plot}
\alias{TAPS}


\title{
Tumor Aberration Prediction Suite
}

\description{
Tumor Aberration Prediction Suite. Version 1.9, April 2013

This package contains all functions used to run (and used by) Tumor Aberration Prediction Suite (TAPS).
(Developed September 2009 --> January 2011, packaged 2013)
Updates: Check http://patchwork.r-forge.r-project.org/ as we intend to extend some functionality.
}


\usage{
TAPS_plot(directory=NULL,autoEstimate=FALSE,bin=400,cores=1,matched=FALSE,allelePeaks=FALSE)
}

\arguments{
	%\item{samples}{
	%Default is 'all'. Alternatively, supply a (boolean, integer or sample name) vector of included samples.
	%}
	\item{directory}{
	Default is NULL. If NULL getwd() is used. Specifying to a specific samples directory will run TAPS_plot on that directory.
	Specifying a directory containing one or more subdirectories that are samples (and not any other subdirectories
	or TAPS_plot will error when trying to run them!) will iteratively run TAPS_plot on all samples.
	}
%	\item{xlim}{
%	Default is c(-1,2).
%	}
	\item{autoEstimate}{
	Default is FALSE. Development stage, do not use.
	}
	\item{bin}{
	Default is 400.
	}
	\item{cores}{
	Default is 1. Only supported on unix based systems. Set amount of threads/cores to be used.
	}
	\item{matched}{
	Default is FALSE. Set to TRUE if a matched normal has been used to remove homozygous snp allele frequencies. 
	}
	\item{allelePeaks}{
	Default is FALSE. Set to TRUE if SNP values are allele peaks, not B-allele frequency. 
	}
}

\details{
Instructions, allele-specific copy number analysis \cr

Note to users:\cr
TAPS is a tool to help you actively investigate genomic aberrations of the most complex tumor samples.
TAPS visualizes samples using log-ratio and allelic imbalance ratio. If the sample is deemed suitable,
 it can be subjected to automated copy number calling.

Data: \cr
Validated automatic copy number calling is available for Affymetrix SNP6, 250k/500k data and cytoscan HD.
Automatic copy number calling is suitable but unvalidated for OncoScan (an Affymetrix SNP array).
If you have some other microarray data please contact us, it is possible that there is a way to make it work.

Data preparation:\cr
Nexus (BioDiscovery): find your "Samples" folder after segmentation (SNP-rank) and normalization.

Other: \cr
Store Log-ratio as probes.txt and allele frequency as snps.txt
in one folder per sample. Do not store other folders than sample folders in your "Samples" folder.
The 'Array' column in these files is uninformative and may be omitted.
Optionally (recommended) also put segments.txt in your sample folder.

Workflow:\cr
1. From the folder containing your samples (sample folders) run TAPS_plot(). \cr
2. Investigate the scatter plots generated in your sample folders, using human interpretation, TAPS_estimates() and TAPS_click(), or just TAPS_click(). \cr
If using human interpretation:
Find and open the file "SampleData.csv". \cr
For each sample, enter an interpretation of Log-Ratio @ copy number 1-3 ("cn1-3", two+ of the three values are required) and the allelic imbalance ratio of LOH and cn2 ("loh"). Save the file. \cr

5. Run TAPS_call(). \cr
6. Inspect the karyotype_check images, and the new chromosome-wise images. \cr
7. If all looks reasonable, you will find good copy number estimates in 'Copynumbers.csv'. \cr
9. Be wary of the result on sex chromosomes which may be difficult to auto-interpret. \cr
10. Watch all images for signs of segmentation failure and tumor cell heterogeneity. \cr
}


\references{
http://www.biomedcentral.com/content/pdf/gb-2011-12-10-r108.pdf
}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Hanna Goransson Kultima, \email{hanna.goransson.kultima@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se} \cr
Bjorn Viklund, \email{bjorn.viklund@medsci.uu.se}
}
