This is a demo or testrun for patchworkCG to help you verify that you have installed
 patchworkCG correctly and that it is working.<br /><br />

It builds upon the knowledge gained from the other tabs, installation to results,
and as such is not overly detailed. Read the other tabs before doing the demo!<br /><br />

 Used in this demo is the 
<a href="ftp://ftp2.completegenomics.com/Cancer_pairs/ASM_Build37_2.0.0/HCC2218/GS00258-DNA_B03/" style="text-decoration:none;"
 target="_blank">publically available cancer data set HCC2218 of CompleteGenomics</a>.
  <br /><br />

We will go through the exact commands applicable in a typical run of this data and
 some plots that you can compare your own results with. <br /><br />

 Once you have downloaded the data and unpacked it you should have a file structure
 resembling the screenshot below. The most important files, which are used in
 patchworkCG, have been highlighted. <br />
 <img src="css/img/file_struct1.png" alt="CompleteGenomics ASM"
title="ASM File structure"> <br /><br />

Now we are ready to start R and issue all the necessary commands, as seen in the
install and execution tabs.

<pre style="height:500px;overflow:scroll;">
#For a working directory I chose just above the ASM folder.

$  GS00258-DNA_B03  ls
ASM                               idMap-HCC2218-H-200-37-ASM-T1.tsv
HCC2218-H-200-37-ASM              version
LIB

#Starting R

$  GS00258-DNA_B03  R

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

#Install package

> install.packages("patchworkCG", repos="http://R-Forge.R-project.org",type="source")
trying URL 'http://R-Forge.R-project.org/src/contrib/patchworkCG_1.0.tar.gz'
Content type 'application/x-gzip' length 669849 bytes (654 Kb)
opened URL
==================================================
downloaded 654 Kb

* installing *source* package 'patchworkCG' ...
** R
** data
** preparing package for lazy loading
** help
*** installing help indices
** building package indices ...
** testing if installed package can be loaded

* DONE (patchworkCG)

The downloaded packages are in
	'/private/var/folders/zm/8v842d8n4172k1yc0589p8sw0000gn/T/RtmprYh5Mq/downloaded_packages'

#Load Library

> library(patchworkCG)

#Execute patchwork.CG.plot

> patchwork.CG.plot(path="ASM/",name="Demo_") 
Reading files from ASM folder 
File input complete 
Performing calculations 
Calculations complete 
Saving objects to CG.Rdata 
Initiating Plotting 
Plotting Complete 
patchwork.CG.plot Complete 
 Please read documentation on running patchwork.CG.copynumbers 

#Execute patchwork.CG.copynumbers (See execution section for aquiring arguments)
#Side note: this is a tricky sample as there is no copy number 2 with loss
#of heterozygosity, but as copynumber 1 has allelic imbalance around 1 and copy number 3
#also has high allelic imbalance we can estimate that had it existed copy number 2
#would have had very high allelic imbalance. So we approximate it to 0.96.

> patchwork.CG.copynumbers(cn2=0.9,delta=0.5,het=0.25,hom=0.96) 
Warning message:
NAs introduced by coercion 

#Quit

> q()
Save workspace image? [y/n/c]: n
$  GS00258-DNA_B03  ls
ASM                               	Demo__KaCh_chr1.png
CG.Rdata                          	Demo__KaCh_chr10.png
copynumbers.Rdata                       Demo__KaCh_chr11.png
copynumbers_KaChCN_chr1.png             Demo__KaCh_chr12.png
copynumbers_KaChCN_chr10.png            Demo__KaCh_chr13.png
copynumbers_KaChCN_chr11.png            Demo__KaCh_chr14.png
copynumbers_KaChCN_chr12.png            Demo__KaCh_chr15.png
copynumbers_KaChCN_chr13.png            Demo__KaCh_chr16.png
copynumbers_KaChCN_chr14.png            Demo__KaCh_chr17.png
copynumbers_KaChCN_chr15.png            Demo__KaCh_chr18.png
copynumbers_KaChCN_chr16.png            Demo__KaCh_chr19.png
copynumbers_KaChCN_chr17.png            Demo__KaCh_chr2.png
copynumbers_KaChCN_chr18.png            Demo__KaCh_chr20.png
copynumbers_KaChCN_chr19.png            Demo__KaCh_chr21.png
copynumbers_KaChCN_chr2.png             Demo__KaCh_chr22.png
copynumbers_KaChCN_chr20.png            Demo__KaCh_chr3.png
copynumbers_KaChCN_chr21.png            Demo__KaCh_chr4.png
copynumbers_KaChCN_chr22.png            Demo__KaCh_chr5.png
copynumbers_KaChCN_chr3.png             Demo__KaCh_chr6.png
copynumbers_KaChCN_chr4.png             Demo__KaCh_chr7.png
copynumbers_KaChCN_chr5.png             Demo__KaCh_chr8.png
copynumbers_KaChCN_chr6.png             Demo__KaCh_chr9.png
copynumbers_KaChCN_chr7.png             Demo__KaCh_chrX.png
copynumbers_KaChCN_chr8.png             Demo__KaCh_chrY.png
copynumbers_KaChCN_chr9.png             Demo__karyotype.png
copynumbers_KaChCN_chrX.png             HCC2218-H-200-37-ASM
copynumbers_KaChCN_chrY.png             LIB
copynumbers_Ka_check.png                idMap-HCC2218-H-200-37-ASM-T1.tsv
Copynumbers.csv                   	version
</pre>

If successfull your working directory should now contain:
<ul class="checks">
	<li>50 plots</li>
	<li>CG.Rdata</li>
	<li>copynumbers.Rdata</li>
	<li>Copynumbers.csv</li>
</ul>
<br />
Here are three of the plots generated during the demo run to validate your plots with. <br /><br />


<a href="css/img/Demo__KaCh_chr1.png" target="_blank" style="text-decoration:none;">
	Demo__KaCh_chr1.png</a><br />
<!-- <img src="css/img/Demo__KaCh_chr1.png" alt="Output patchwork.CG.plot"
title="Chromosome 1, patchwork.CG.plot" width=880px> -->


<a href="css/img/copynumbers_KaChCN_chr1.png" target="_blank" style="text-decoration:none;">
	copynumbers_KaChCN_chr1.png</a><br />

<!-- <img src="css/img/copynumbers_KaChCN_chr1.png" alt="Output1 patchwork.CG.copynumbers"
title="Chromsome 1, patchwork.CG.copynumbers" width=880px> -->

<a href="css/img/copynumbers_Ka_check.png" target="_blank" style="text-decoration:none;">
	copynumbers_Ka_check.png</a><br />

<!-- <img src="css/img/copynumbers_Ka_check.png" alt="Output2 patchwork.CG.copynumbers"
title="Whole Genome, Check" width=880px><br /><br /> -->

<br /> <br /> Please feel free to contact us with any questions, suggestions or feedback!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>




















