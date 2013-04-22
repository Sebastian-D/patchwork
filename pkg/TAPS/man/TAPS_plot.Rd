\name{TAPS_plot}
\alias{TAPS_plot}
\alias{Taps_plot}
\alias{TAPS.plot}
\alias{Taps.plot}
\alias{TAPS}


\title{
Tumor Abberation Prediction Suite
}

\description{
Tumor Aberration Prediction Suite. Version 1.9, April 2013

This package contains all functions used to run (and used by) Tumor Aberration Prediction Suite (TAPS).
Author: Markus Rasmussen @ Uppsala University. (Developed September 2009 --> January 2011, packaged 2013)
Contact: markus.rasmussen@medsci.uu.se or arrayplatform@medsci.uu.se
Updates: Check http://patchwork.r-forge.r-project.org/ as we intend to extend some functionality.
}

\usage{
TAPS_plot(directory=NULL,bin=400)
}

\arguments{
	\item{directory}{
	Default is getwd().
	}
%	\item{xlim}{
%	Default is c(-1,2).
%	}
%	\item{ylim}{
%	Default is c(0,1).
%	}
	\item{bin}{
	Default is 400.
	}
}

\details{
-----Instructions, allele-specific copy number analysis:-------------------------------- \\

Note to users:
TAPS is a tool to help you actively investigate genomic aberrations of the most complex tumor samples.
TAPS visualizes samples using log-ratio and allelic imbalance ratio. If the sample is deemed suitable,
 it can be subjected to automated copy number calling.
Data: Validated automatic copy number calling is available for Affymetrix SNP6 and 250k/500k data.
  Automatic copy number calling is suitable but unvalidated for OncoScan (an Affymetrix SNP array)

Data preparation.
Nexus (BioDiscovery) 5.0: find your "Samples" folder after segmentation (SNP-rank) and normalization.

Other: See example package folders. Store Log-ratio as probes.txt and allele frequency as snps.txt
in one folder per sample. Do not store other folders than sample folders in your "Samples" folder.
The 'Array' column in these files is uninformative and may be omitted.
Optionally (recommended) also put segments.txt in your sample folder. Please inspect the examples
for table structure and headers.

Workflow.
1. From the folder containing your samples (sample folders) run TAPS_plot().
2. Investigate the scatter plots generated in your sample folders
3. To proceed with copy number calls, find and open the file "SampleData.txt"
4. For each sample, enter an interpretation of Log-Ratio @ copy number 2 ("cn2"), the difference in
   Log-Ratio to a deletion ("delta") and the allelic imbalance ratio of CNNLOH ("loh"). Save the file.
5. Run TAPS_call().
6. Inspect the karyotype_check images, and the new chromosome-wise images.
7. If all looks reasonable, you will find good copy number estimates in 'Copynumbers.csv'.
9. Be wary of the result on sex chromosomes which may be difficult to auto-interpret.
10. Watch all images for signs of segmentation failure and tumor cell heterogeneity.
}


\references{
http://www.biomedcentral.com/content/pdf/gb-2011-12-10-r108.pdf
}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se} \cr
Hanna Goransson Kultima, \email{hanna.goransson.kultima@medsci.uu.se} \cr
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
