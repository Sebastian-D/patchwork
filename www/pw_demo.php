To help you verify that you have installed patchwork correctly and that
it is working we will use these files for a testrun. <br /><br />

<a href="http://130.238.204.28/open/patchwork/tumor.bam" style="text-decoration:none;" target="_blank">
tumor.bam (72GB)</a><br />
<a href="http://130.238.204.28/open/patchwork/tumor_pileup" style="text-decoration:none;" target="_blank">
tumor_pileup (441 MB)</a><br />
<a href="http://130.238.204.28/open/patchwork/tumor.bam.bai" style="text-decoration:none;" target="_blank">
tumor.bam.bai (8.7 MB)</a><br />
<a href="http://130.238.204.28/open/patchwork/references/datasolexa.RData" style="text-decoration:none;" target="_blank">
Solexa/Illumina reference (64.8 MB)</a><br /><br />


  We will go through the exact commands applicable in a typical run of this data and
 some plots that you can compare your own results with. <br /><br />

Once you have downloaded the files we start R and issue all the necessary commands.
In this demo the working directory is the same folder where the example files are
located.

<pre style="height:500px;overflow:scroll;">
#My working directory, with the example files

➜  pw_demo  ls
tumor.bam.bai 	  tumor_pileup  
datasolexa.RData  tumor.bam

#Starting R

➜  pw_demo  R

R version 2.14.1 (2011-12-22)
Copyright (C) 2011 The R Foundation for Statistical Computing
ISBN 3-900051-07-0
Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

#Installing DNAcopy

> source("http://bioconductor.org/biocLite.R")
BiocInstaller version 1.2.1, ?biocLite for help
>     biocLite("DNAcopy")
BioC_mirror: 'http://www.bioconductor.org'
Using R version 2.14, BiocInstaller version 1.2.1.
Installing package(s) 'DNAcopy'
trying URL 'http://www.bioconductor.org/packages/2.9/bioc/bin/macosx/leopard/contrib/2.14/DNAcopy_1.28.0.tgz'
Content type 'application/x-gzip' length 461156 bytes (450 Kb)
opened URL
==================================================
downloaded 450 Kb


The downloaded packages are in
	/var/folders/zm/8v842d8n4172k1yc0589p8sw0000gn/T//RtmpVqklVr/downloaded_packages
> 

#Installing patchwork/patchworkData

> install.packages("patchwork", repos="http://R-Forge.R-project.org")
trying URL 'http://R-Forge.R-project.org/bin/macosx/leopard/contrib/2.14/patchwork_2.1.tgz'
Content type 'application/x-gzip' length 933161 bytes (911 Kb)
opened URL
==================================================
downloaded 911 Kb


The downloaded packages are in
	/var/folders/zm/8v842d8n4172k1yc0589p8sw0000gn/T//Rtmpa2Vj4Z/downloaded_packages
> install.packages("patchworkData", repos="http://R-Forge.R-project.org")
trying URL 'http://R-Forge.R-project.org/bin/macosx/leopard/contrib/2.14/patchworkData_2.0.tgz'
Content type 'application/x-gzip' length 91933481 bytes (87.7 Mb)
opened URL
=================================================
downloaded 87.7 Mb


The downloaded packages are in
	/var/folders/zm/8v842d8n4172k1yc0589p8sw0000gn/T//Rtmpa2Vj4Z/downloaded_packages

#Loading the patchwork Libraries

> library(patchwork)

> library(patchworkData)

#Execute patchwork.plot

> patchwork.plot(BamFile="tumor.bam",Pileup="tumor_pileup",reference="datasolexa.RData")
Initiating Allele Data Generation 
Allele Data Generation Complete 
Initiating Read Chromosomal Coverage 
Reading chr1 
Reading chr2 
Reading chr3 
Reading chr4 
Reading chr5 
Reading chr6 
Reading chr7 
Reading chr8 
Reading chr9 
Reading chrX 
Reading chrY 
Reading chr10 
Reading chr11 
Reading chr12 
Reading chr13 
Reading chr14 
Reading chr15 
Reading chr16 
Reading chr17 
Reading chr18 
Reading chr19 
Reading chr20 
Reading chr21 
Reading chr22 
Read Chromosomal Coverage Complete 
Initiating GC Content Normalization 
GC Content Normalization Complete 
Initiating Smoothing 
Smoothing Chromosome: chr1 
Smoothing Chromosome: chr2 
Smoothing Chromosome: chr3 
Smoothing Chromosome: chr4 
Smoothing Chromosome: chr5 
Smoothing Chromosome: chr6 
Smoothing Chromosome: chr7 
Smoothing Chromosome: chr8 
Smoothing Chromosome: chr9 
Smoothing Chromosome: chrX 
Smoothing Chromosome: chrY 
Smoothing Chromosome: chr10 
Smoothing Chromosome: chr11 
Smoothing Chromosome: chr12 
Smoothing Chromosome: chr13 
Smoothing Chromosome: chr14 
Smoothing Chromosome: chr15 
Smoothing Chromosome: chr16 
Smoothing Chromosome: chr17 
Smoothing Chromosome: chr18 
Smoothing Chromosome: chr19 
Smoothing Chromosome: chr20 
Smoothing Chromosome: chr21 
Smoothing Chromosome: chr22 
Smoothing Complete 

**************************************************************************
   The plan to change the data format for CNA object has been postponed   
 in order to ensure backward compatibility with older versions of DNAcopy 
**************************************************************************

Initiating Segmentation 
Note: If segmentation fails to initiate the probable reason is that you have not installed the R package DNAcopy. See patchworks README for installation instructions. 
Analyzing: chr1.p 
Analyzing: chr1.q 
Analyzing: chr2.p 
Analyzing: chr2.q 
Analyzing: chr3.p 
Analyzing: chr3.q 
Analyzing: chr4.p 
Analyzing: chr4.q 
Analyzing: chr5.p 
Analyzing: chr5.q 
Analyzing: chr6.p 
Analyzing: chr6.q 
Analyzing: chr7.p 
Analyzing: chr7.q 
Analyzing: chr8.p 
Analyzing: chr8.q 
Analyzing: chr9.p 
Analyzing: chr9.q 
Analyzing: chrX.p 
Analyzing: chrX.q 
Analyzing: chrY.p 
Analyzing: chrY.q 
Analyzing: chr10.p 
Analyzing: chr10.q 
Analyzing: chr11.p 
Analyzing: chr11.q 
Analyzing: chr12.p 
Analyzing: chr12.q 
Analyzing: chr13.q 
Analyzing: chr14.q 
Analyzing: chr15.q 
Analyzing: chr16.p 
Analyzing: chr16.q 
Analyzing: chr17.p 
Analyzing: chr17.q 
Analyzing: chr18.p 
Analyzing: chr18.q 
Analyzing: chr19.p 
Analyzing: chr19.q 
Analyzing: chr20.p 
Analyzing: chr20.q 
Analyzing: chr21.p 
Analyzing: chr21.q 
Analyzing: chr22.q 
Segmentation Complete 
Initiating Segment data extraction (Medians and AI) 
Segment data extraction Complete 
 
 
Saving information objects needed for patchwork.copynumbers in tumor_copynumbers.Rdata 
 
 
Initiating Plotting 
Plotting Complete 
Shutting down..... 
Warning messages:
1: In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'tumor_pile.alleles.Rdata', probable reason 'No such file or directory'
2: In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'tumor_data.Rdata', probable reason 'No such file or directory'
3: In if (system(paste("python ", getwd(), "/.python/Check.py", sep = ""),  :
  the condition has length > 1 and only the first element will be used
4: In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'tumor_smoothed.Rdata', probable reason 'No such file or directory'
5: In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'tumor_Segments.Rdata', probable reason 'No such file or directory'

#Note that these warning messages are completely fine. See "execution tab" for more info.

#Execute patchwork.copynumbers
#(For guide on how to discern arguments see "execution tab")

> patchwork.copynumbers(cn2=0.68,delta=0.23,het=0.42,hom=0.78)

#Quit

> q()
Save workspace image? [y/n/c]: n

➜  pw_demo ls
copynumbers__karyotype_check.png    tumor_Segments.Rdata
copynumbers__karyotypeCN.chr10.png  tumor.bam.karyotype.png
copynumbers__karyotypeCN.chr11.png  tumor_pileup
copynumbers__karyotypeCN.chr12.png  tumor_smoothed.Rdata
copynumbers__karyotypeCN.chr13.png  tumor.bam
copynumbers__karyotypeCN.chr14.png  tumor.bam.bai
copynumbers__karyotypeCN.chr15.png  tumor.bam_karyotype.chr10.png
copynumbers__karyotypeCN.chr16.png  tumor.bam_karyotype.chr11.png
copynumbers__karyotypeCN.chr17.png  tumor.bam_karyotype.chr12.png
copynumbers__karyotypeCN.chr18.png  tumor.bam_karyotype.chr13.png
copynumbers__karyotypeCN.chr19.png  tumor.bam_karyotype.chr14.png
copynumbers__karyotypeCN.chr1.png   tumor.bam_karyotype.chr15.png
copynumbers__karyotypeCN.chr20.png  tumor.bam_karyotype.chr16.png
copynumbers__karyotypeCN.chr21.png  tumor.bam_karyotype.chr17.png
copynumbers__karyotypeCN.chr22.png  tumor.bam_karyotype.chr18.png
copynumbers__karyotypeCN.chr2.png   tumor.bam_karyotype.chr19.png
copynumbers__karyotypeCN.chr3.png   tumor.bam_karyotype.chr1.png
copynumbers__karyotypeCN.chr4.png   tumor.bam_karyotype.chr20.png
copynumbers__karyotypeCN.chr5.png   tumor.bam_karyotype.chr21.png
copynumbers__karyotypeCN.chr6.png   tumor.bam_karyotype.chr22.png
copynumbers__karyotypeCN.chr7.png   tumor.bam_karyotype.chr2.png
copynumbers__karyotypeCN.chr8.png   tumor.bam_karyotype.chr3.png
copynumbers__karyotypeCN.chr9.png   tumor.bam_karyotype.chr4.png
copynumbers__karyotypeCN.chrX.png   tumor.bam_karyotype.chr5.png
copynumbers__karyotypeCN.chrY.png   tumor.bam_karyotype.chr6.png
Copynumbers.csv 		    tumor.bam_karyotype.chr7.png
tumor_copynumbers.Rdata                   tumor.bam_karyotype.chr8.png
tumor_data.Rdata                          tumor.bam_karyotype.chr9.png
tumor_datasolexa.RData                    tumor.bam_karyotype.chrX.png
pile.alleles                        tumor.bam_karyotype.chrY.png
tumor_pile.alleles.Rdata                  
</pre>

If successfull your working directory should now contain:
<ul class="checks">
	<li>50 plots</li>
	<li>tumor_copynumbers.Rdata</li>
	<li>tumor_data.Rdata</li>
	<li>pile.alleles</li>
	<li>tumor_pile.alleles.Rdata</li>
	<li>tumor_Segments.Rdata</li>
	<li>tumor_smoothed.Rdata</li>
</ul>
<br />

Note, as explained in the execution tab, that the files are named "tumor_" because the input .bam is named
tumor.bam.
<br /><br />


Here are three of the plots to validate your plots with.
<br /><br />

<a href="css/img/tumor.bam_karyotype.chr18.png" target="_blank" style="text-decoration:none;">
	tumor.bam_karyotype.chr18.png</a><br />

<a href="css/img/copynumbers__karyotypeCN.chr18.png" target="_blank" style="text-decoration:none;">
	copynumbers__karyotypeCN.chr18.png</a><br />


<a href="css/img/copynumbers__karyotype_check.png" target="_blank" style="text-decoration:none;">
	copynumbers__karyotype_check.png</a><br />

<br /> <br /> Please feel free to contact us with any queries!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>


