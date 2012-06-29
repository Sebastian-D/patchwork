<!-- Flowchart -->
<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowchartCG.png" 
alt="Flowchart of processes and data for patchworkCG" title="Flowchart of data and processes for patchworkCG">
</div><br />

<h4>patchwork.CG.plot()</h4>

patchwork.CG.plot() is a function that visualizes allele-specific copy numbers from the
whole-genome sequenced data of tumors. <br /><br />
It is recommended that you run patchworkCG from a "clean" working directory.
This prevents the risk of having files write over eachother when
 you run multiple different samples.<br /><br />

In the R environment load the patchworkCG library: <br />

<pre>
	library(patchworkCG)
</pre>

Read the excellent documentation for patchwork.CG.plot(): <br />

<pre>
	?patchwork.CG.plot
</pre>

Excerpt from ?patchwork.CG.plot regarding usage and arguments: <br />

<pre>
	Usage:

     	patchwork.CG.plot(path,name='CG_sample',manual_file_input=FALSE,masterVarBeta=NULL,somaticCnvSegments=NULL,
     			depthOfCoverage=NULL)
     
	Arguments:

	    path: Path to the ASM folder. One of the subfolders of your
	          completegenomics directory.
	    name: Default is 'CG_sample'. The name you wish associated with
	          the plots that will be generated.
	manual_file_input: Default is FALSE. If you set it to TRUE you will be
	          prompted to provide path and filename for the masterVarBeta,
	          somaticCnvSegmentsNondiploid and depthOfCoverage file.
	masterVarBeta: Path to and COMPLETE NAME of the masterVarBeta file,
	          should you wish to implement it directly. Useful if you want
	          to run the program on multiple samples and wish to create a
	          script going from file to file. This saves you from having to
	          conserve the CompleteGenomics file structure.  Default is
	          NULL.
	somaticCnvSegments: Path to and COMPLETE NAME of the somaticCnvSegments
	          file, should you wish to implement it directly. Useful if you
	          want to run the program on multiple samples and wish to
	          create a script going from file to file. This saves you from
	          having to conserve the CompleteGenomics file structure.
	          Default is NULL.
	depthOfCoverage: Path to and COMPLETE NAME of the depthOfCoverage file,
	          should you wish to implement it directly. Useful if you want
	          to run the program on multiple samples and wish to create a
	          script going from file to file. This saves you from having to
	          conserve the CompleteGenomics file structure.  Default is
	          NULL.
</pre>

Example run of patchwork.CG.plot() where output has been left for display: <br />

<pre>
	patchwork.CG.plot(path="path/to/ASM/")

		Reading files from ASM folder 
		File input complete 
		Performing calculations 
		Calculations complete 
		Saving objects to CG.Rdata 
		Initiating Plotting 
		Plotting Complete 
		patchwork.CG.plot Complete 
		 Please read documentation on running patchwork.CG.copynumbers 
</pre>

Once it is complete there should be several plots in your working directory as
well as the file "CG.Rdata". <br /><br />

Execution should take no longer than half an hour, depending on your system.<br /><br />

<h4>patchwork.CG.copynumbers()</h4>

patchwork.CG.copynumbers() is a function that determines which relationship between coverage
  and allelic imbalance signifies which copy number and allele ratio for each segment.<br /><br />

Running patchwork.CG.copynumbers() requires the aforementioned CG.Rdata file, the function
will find it automatically if it is in your working directory. There is also
the option to point to the file using the "CGfile" argument of the function.<br /><br />

Read the documentation for patchwork.CG.copynumbers(): <br />

<pre>
	?patchwork.CG.copynumbers
</pre>

Excerpt from ?patchwork.CG.copynumbers regarding usage and arguments: <br />

<pre>
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
</pre>

To infer the arguments for patchwork.CG.copynumbers() you will need to look at one of the chromosomal
 plots generated using patchwork.CG.plot(). Here is the topmost section of one such plot with some added bars
  and text for tutorial purposes. The plot shows the whole genome in grey and the chosen chromosome, not
  important in this guide, in red to blue gradient. From this whole genome picture you can then use the
  arrangement of the clustered areas to estimate copy numbers and allele-specific information. <br /><br />

<img id="fig1cg" src="css/img/Fig1pw_tut.png" 
alt="Input arguments patchwork.CG.copynumbers" title="Figure 1 patchworkCG"
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

There is more information on interpreting the plots in the results tab which
may aid you in determining the parameters. <br /><br />

An example run of patchwork.CG.copynumbers(): <br />

<pre>
	patchwork.CG.copynumbers(cn2=0.8, delta=0.28, het=0.21, hom=0.79, CGfile="path/to/CG.Rdata")
</pre>

This will generate an additional set of plots in your working directory. <br /><br />


<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>




























