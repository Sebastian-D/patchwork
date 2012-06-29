<!-- Flowchart -->
<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowchartPWexec.png" 
alt="Flowchart of processes and data for patchwork" title="Flowchart of data and processes for patchwork">
</div><br />

<h4>patchwork.plot()</h4>

patchwork.plot() is a function that visualizes allele-specific copy numbers from the
 whole-genome sequenced data of tumors.<br /><br />

It is recommended that you run patchwork from a "clean" working directory to avoid the risk
 of having files write over eachother if you run multiple different samples.<br /><br />

Execution may take quite a while depending on the size of your sample!
So if possible run it on a dedicated computer.<br /><br />

Initiate the R environment and load the patchwork and patchworkData libraries: <br />

<pre>
	library(patchwork)
	library(patchworkData)
</pre>

Read the excellent documentation for patchwork.plot(): <br />

<pre>
	?patchwork.plot
</pre>

Excerpt from ?patchwork.plot regarding usage and arguments: <br />

<pre>
	Usage:

	     patchwork.plot(BamFile,Pileup,Reference=NULL,Normal.bam=NULL,Normal.pileup=NULL,Alpha=0.0001,SD=1)
	     
	Arguments:

	 BamFile: Aligned, Sorted Bam-file.
	  Pileup: Pileup file generated through samtools -vcf reference.fasta
	          bamfile > outfile.
	Reference: Default is NULL.  Path to a reference file that can be
	          created using patchwork.createreference() or downloaded from
	          patchworks homepage.
	Normal.bam: Default is NULL.  The matched normal sample of the your
	          BamFile.
	Normal.pileup: Default is NULL.  The pileup of your normal sample.
	   Alpha: Default 0.0001, change if you want to try to get a better
	          segmentation from patchwork.segment().  From DNAcopy
	          (?segment): alpha: significance levels for the test to accept
	          change-points.
	      SD: Default 1, change if you want to try to get a better
	          segmentation from patchwork.segment().  From DNAcopy
	          (?segment): undo.SD: the number of SDs between means to keep
	          a split if undo.splits="sdundo".
</pre>


Example run of patchwork.plot() where output has been left for display: <br />

<pre style="height:200px;overflow:scroll;">
	patchwork.plot(BamFile="patchwork.example.bam",Pileup="patchwork.example.pileup",Reference="../HCC1954/datasolexa.RData")
		Initiating Allele Data Generation
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
		 
		 
		Saving information objects needed for patchwork.copynumbers in (prefix)_copynumbers.Rdata 
		 
		 
		Initiating Plotting 
		Plotting Complete 
		Shutting down..... 
		Warning messages:
		(Some warning messages about trying to load the files in the checklist below is not uncommon)
</pre>

Your working directory should now have the plots generated from the function, 1 overview
 plot and 24 chromosomal plots. <br /><br />

The working directory should also contain these files, where (prefix) is the name of your BAM file
 without the .bam file ending : <br />
<ul class="checks">
	<li>(prefix)_copynumbers.Rdata</li>
	<li>(prefix)_data.Rdata</li>
	<li>(prefix)_pile.alleles.Rdata</li>
	<li>(prefix)_Segments.Rdata</li>
	<li>(prefix)_smoothed.Rdata</li>
	<li>pile.alleles</li>
</ul>

These were created for swifter re-runs of the function should something unforseen happen during execution.
For example if something goes wrong during a final step it would be a terrible hassle to read all
the chromosomes again. <br /><br />

<div style="text-align: center;">
<b>
TO COMPLETELY RE-RUN PATCHWORK FROM SCRATCH ALL OF THE ABOVE FILES NEED TO BE DELETED FROM YOUR WORKING DIRECTORY
</b>
</div>
<br />

 <h4>patchwork.copynumbers()</h4>


 patchwork.copynumbers() is a function that determines which relationship between coverage
  and allelic imbalance signifies what copy number and allele ratio for each segment. <br /><br />

The only file you absolutely must have in your working for the next part of execution is (prefix)_copynumbers.Rdata. <br /><br />

Read the outstanding documentation for patchwork.copynumbers():<br />

<pre>
	?patchwork.copynumbers
</pre>

Excerpt from ?patchwork.copynumbers regarding usage and arguments: <br />

<pre>
	Usage:

	     patchwork.copynumbers(CNfile,cn2,delta,het,hom,maxCn=8,ceiling=1,forcedelta=F)
	     
	Arguments:

	    CNfile: The name and path of your copynumbers file, generated from patchwork.plot().
	    		Example Myfile_copynumbers.Rdata.
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

</pre>


To infer the arguments for patchwork.copynumbers() you will need to look at one of the chromosomal
 plots generated using patchwork.plot(). Here is the topmost section of one such plot with some added bars
  and text for tutorial purposes. The plot shows the whole genome in grey and the chosen chromosome, not
  important in this guide, in red to blue gradient. From this whole genome picture you can then use the
  arrangement of the clustered areas to estimate copy numbers and allele-specific information. <br /><br />

<img id="fig1pw" src="css/img/Fig1pw_tut.png" alt="Input arguments patchwork.copynumbers" title="Figure 1 patchwork"
 style="width:880px"><br /><br />


The <b>cn2</b> argument is the position of copy number 2. In this example cn2 is ~0.8.<br />
The <b>delta</b> argument is the difference between two copynumbers on the coverage axis. In this example we
used copy number 2 and copy number 3, but we could as well have used copy number 2 and the unmarked copy number
1 to the left of copy number 2. Delta does not take negative argument input. In this example delta is ~0.28.<br />
The <b>het</b> argument is the position of heterozygous copy number 2 on the allelic imbalance axis. In this example
het is ~0.21.<br />
The <b>hom</b> argument is the position of homozygous, loss-of-heterozygosity, copy number 2 on the allelic imbalance
axis. In this example hom is ~0.79. <br /><br />

Usually you do not need to bother with <b>maxCn</b> and <b>ceiling</b> arguments. <br /><br />

An example run of patchwork.copynumbers(): <br />

<pre>
	patchwork.copynumbers(cn2=0.8, delta=0.28, het=0.21, hom=0.79)
</pre>

This will generate an additional set of plots in your working directory. <br /><br />

In the results tab we will look at interpreting the data.

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>

























































