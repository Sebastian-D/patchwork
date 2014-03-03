<!-- Flowchart -->
<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowchartTAPS.jpg" 
alt="Flowchart of processes and data for patchwork" title="Flowchart of data and processes for patchwork">
</div><br />

<h4>TAPS_plot()</h4>

TAPS_plot() is a function that visualizes allele-specific copy numbers from microarray
data of tumors. It can either be run in a directory containing multiple sample folders or
pointing to a specific sample, which would quite simply be a folder containing no subfolders
and ideally the files needed to execute.
<br /><br />

The execution of TAPS_plot() should be pretty quick for samples that have been processed
by Nexus, as the segmentation step is already completed. For other samples it may take
about an hour. <br /><br />

To run TAPS_plot(), start by initiating the R environment and loading the TAPS library: <br />

<pre>
	library(TAPS)
</pre>

Read the documentation for TAPS_plot(): <br />

<pre>
	?TAPS_plot
</pre>

Excerpt from ?TAPS_plot regarding usage and arguments: <br />

<pre>
Usage:

     TAPS_plot(directory=NULL,bin=400)
     
Arguments:

directory: Default is getwd(). Specifying to a specific samples
          directory will run TAPS_plot on that directory.  Specifying a
          directory containing one or more subdirectories that are
          samples (and not any other subdirectories or TAPS_plot will
          error when trying to run them!) will iteratively run
          TAPS_plot on all samples.

     bin: Default is 400.
</pre>

The only thing TAPS_plot() needs is the directory argument. If, as mentioned above, you have
a folder containing several sample folders then specifying "directory" to that main folder
will execute all the sample subfolders.

<img id="main" src="css/img/Main_folder.png" alt="Main folder" title="Main folder input"
 style="border:1px solid black"><br /><br />

Specifying "directory" to a single sample folder will execute only that sample. For example
the folders Sample1 and Sample2.
<br /><br />
<b>Nexus</b><br />
 <img id="Nexus" src="css/img/Nexus.png" alt="Nexus" title="Nexus input"
 style="border:1px solid black">
 <br />
<b>ChAS</b><br />
<img id="Chas" src="css/img/Chas.png" alt="Chas" title="Chas input"
style="border:1px solid black">
<br /><br />

Note that specifying any folder containing subfolders will interpret the subfolders as the sample folders.
<br /><br />

Here is an example execution of TAPS_plot() using the current working directory, ".", as input
and all the output left for display. The 6295 sample was processed from Nexus. TAPS_plot()
processed sample 6295 completely then went into the next folder, "NonSampleFolder", which is not
a sample folder and thus threw a bunch of warning messages and failed. <br /><br />

<pre style="height:200px;overflow:scroll;">
TAPS_plot(directory=".")
   ..loading 6295 ..found probes.txt ..found snps.txt ..found segments.txt ..processing..plotting.
  ..done
   ..loading NonSampleFolder ..processingUsing DNAcopy to create segments:
  Warning messages:
  1: In file(file, "rt") :
    cannot open file 'probes.txt': No such file or directory
  2: In file(file, "rt") :
    cannot open file '_probes.txt': No such file or directory
  3: In file(file, "rt") :
    cannot open file 'snps.txt': No such file or directory
  4: In file(file, "rt") :
    cannot open file '_snps.txt': No such file or directory
  5: In is.na(Log2$Value) :
    is.na() applied to non-(list or vector) of type 'NULL'
  6: In is.na(alf$Value) :
    is.na() applied to non-(list or vector) of type 'NULL'
  7: In file(file, "rt") :
    cannot open file 'segments.txt': No such file or directory
  8: In file(file, "rt") :
    cannot open file '_segments.txt': No such file or directory
</pre>

This will have generated one "&lt;samplename&gt;_overview.jpeg" in the main directory (should
there be one), several plots and .Rdata files in the sample folder and a file named
simply "SampleData.csv" which will be covered shortly. <br /><br />

<hr class="alt1" />

<h4>User Evaluation</h4>

To infer the arguments for SampleData.csv that TAPS_call() uses you will need to look at one of
the chromosomal plots generated using TAPS_plot(). The structure and relationships in the plot
can be interpreted to figure out the most probable locations of the allele-specific copynumbers.
<br /><br />

For information to help you understand the structure of the plot and how the axis allelic imbalance
 and average log-ratio can be interpreted,
click <a href="AI_ALR.php" target="_blank" style="text-decoration:none;">here</a>.
<br /><br />

Now lets take a look at the structure and placement of clusters on the whole genome plot.<br /><br />
As average log-ratio increases, copy numbers increases. As Allelic imbalance increases, heterozygosity decreases,
that is to say that clusters with high allelic imbalance will be homozygous.<br /><br />
What do we expect a hypothetical plots arrangement of clusters to look like?
The average ploidy of the sample will be 0 on average log-ratio axis as it is normalized.
The sample may be highly rearranged but quite often this is a starting point for finding
copynumber 2 or copynumber 3. As there will be less reads covering copynumber 1 in the sample
than higher copynumbers and copynumber 1 cannot have different allele constitutions, by
its very nature of being one allele, copynumber 1 will be represented by a single cluster far to the left
on average log-ratio axis when compared to the other clusters.<br /><br />
So what if we do not have any copynumber 1 in the sample? Then perhaps the far left
of the plot will be occupied by two clusters, indicating the homozygous, LOH, and heterozygous states
of copynumber 2.
It then stands to reason that the next cluster we will encounter, again; moving from left to right
on average log-ratio axis, will be copynumber 3.<br />
In the same way we would expect copynumber 2 to follow
copynumber 1 in the previous scenario. Note also that the lowest allelic imbalance of copynumber 3
MUST be higher than the lowest allelic imbalance of copynumber 2 as copynumber 2 can have one of each allele
while copynumber 3 must, at lowest allelic imbalance, have two of one allele and one of the other.<br />

The X chromosome is shown as "X"'s instead of gradient circles and may also be far left, do not let it
disturb your assessment of copynumber levels as it can be misleading. <br /><br />

It is with reasoning such as this, looking at the plot and how the clusters are arranged and what cluster constitutions are
physically possible, that we can determine the allele-specific copy numbers in our sample.<br /><br />

Here is the topmost section of one such plot with some added bars and text for tutorial purposes. The plot
shows the whole sample in grey and the chosen chromosome, not important in this guide, in red to blue gradient.
From this whole sample picture you can then use the arrangement of the clustered areas to estimate copy numbers
and allele-specific information. <br /><br />


<img id="fig1taps" src="css/img/TAPS_call_help2.png" alt="Input arguments TAPS_call" title="Figure 1 TAPS"
 style="width:880px"><br /><br />

 The <b>cn1</b>, <b>cn2</b> and <b>cn3</b> arguments are the positions of copy number 1, 2 and 3 on
the "Average log-ratio" axis. In this example cn2 is ~ -0.15. You need only provide any two of these values. <br />
The <b>loh</b> argument is the position of homozygous copy number 2 (LOH) on the "Allelic imbalance" axis.
In this example loh is ~ 0.75.
<br /><br />

We have developed two functions, TAPS_estimate() and TAPS_click, that you can use if you do not want to fill SampleData.csv by hand.
This is especially useful if you are running multiple samples. We would recommend that you always use TAPS_click().
However, you will still have to be able to interpret the plots!
<br /><br />

<h3>TAPS_estimate()</h3>

TAPS_estimate() attempts to estimate the SampleData.csv arguments for all samples in the current working directory.
This information is hard to automate therefore it will not always do so correctly, but it will still simplify the process
when running TAPS_click() where you can either choose to accept TAPS_estimates() evaluation or redo it yourself.
<br /><br />

<h3>TAPS_click()</h3>

<!-- Script for a slideshow -->
<SCRIPT LANGUAGE="JavaScript">

var num=1
img1 = new Image ()
img1.src = "css/img/TAPS_call_help2.png"
img2 = new Image ()
img2.src = "css/img/70_HCC1187__KaCh_chr1.png"
img3 = new Image ()
img3.src = "css/img/70_HCC1187_check.png"
img4 = new Image ()
img4.src = "css/img/Allelic_plot.png" 

function slideshowNext()
{
num=num+1
if (num==5)
{num=1}
document.s1.src=eval("img"+num+".src")
document.s2.value=eval("text"+num)
}

function slideshowBack()
{
num=num-1
if (num==0)
{num=4}
document.s1.src=eval("img"+num+".src")
document.s2.value=eval("text"+num)
}

</SCRIPT>

<CENTER>
<IMG SRC="css/img/TAPS_call_help2.png" NAME="mypic" BORDER=0 HEIGHT="100%" WIDTH="880px">
<p>

<A HREF="JavaScript:slideshowBack()" style="text-decoration:none;"> Back</A>

<A HREF="JavaScript:slideshowNext()" style="text-decoration:none;"> Next</A> 




<hr class="alt1" />

<h4>TAPS_call()</h4>

The function TAPS_call() uses the relationship between log-ratio
and allelic imbalance to assign copy number and allele ratio for each segment.
<br /><br />

TAPS_plot() must have been executed before TAPS_call() can be. This is because the .Rdata files
generated during TAPS_plot() will be automatically read by TAPS_call() and you need to evaluate the plots
generated by TAPS_plot() to be able to give TAPS_call() meaningful input.  
<br /><br />

TAPS_call() also has "directory" as its main argument, the parameters are read from the file
"SampleData.csv". If you desire to execute TAPS_call() on many samples, add the arguments to the
main directory SampleData.csv. If you want to execute TAPS_call() on just one sample, add the arguments
for that sample to SampleData.csv and specify that sample using the "directory" argument of TAPS_call().
<br /><br />

An example run of TAPS_call(), where "mainsamplefolder" contains SampleData.csv: <br />

<pre>
	TAPS_call(directory="path/to/mainsamplefolder")
</pre>

<!-- To infer the arguments for SampleData.csv that TAPS_call() uses you will need to look at one of
the chromosomal plots generated using TAPS_plot(). The structure and relationships in the plot
can be interpreted to figure out the most probable locations of the allele-specific copynumbers.
<br /><br />

For information to help you understand the structure of the plot and how the axis allelic imbalance
 and average log-ratio can be interpreted,
click <a href="AI_ALR.php" target="_blank" style="text-decoration:none;">here</a>.
<br /><br />

Now lets take a look at the structure and placement of clusters on the whole genome plot.<br /><br />
As average log-ratio increases, copy numbers increases. As Allelic imbalance increases, heterozygosity decreases,
that is to say that clusters with high allelic imbalance will be homozygous.<br /><br />
What do we expect a hypothetical plots arrangement of clusters to look like?
The average ploidy of the sample will be 0 on average log-ratio axis as it is normalized.
The sample may be highly rearranged but quite often this is a starting point for finding
copynumber 2 or copynumber 3. As there will be less reads covering copynumber 1 in the sample
than higher copynumbers and copynumber 1 cannot have different allele constitutions, by
its very nature of being one allele, copynumber 1 will be represented by a single cluster far to the left
on average log-ratio axis when compared to the other clusters.<br /><br />
So what if we do not have any copynumber 1 in the sample? Then perhaps the far left
of the plot will be occupied by two clusters, indicating the homozygous, LOH, and heterozygous states
of copynumber 2.
It then stands to reason that the next cluster we will encounter, again; moving from left to right
on average log-ratio axis, will be copynumber 3.<br />
In the same way we would expect copynumber 2 to follow
copynumber 1 in the previous scenario. Note also that the lowest allelic imbalance of copynumber 3
MUST be higher than the lowest allelic imbalance of copynumber 2 as copynumber 2 can have one of each allele
while copynumber 3 must, at lowest allelic imbalance, have two of one allele and one of the other.<br />

The X chromosome is shown as "X"'s instead of gradient circles and may also be far left, do not let it
disturb your assessment of copynumber levels as it can be misleading. <br /><br />

It is with reasoning such as this, looking at the plot and how the clusters are arranged and what cluster constitutions are
physically possible, that we can determine the allele-specific copy numbers in our sample.<br /><br />

Here is the topmost section of one such plot with some added bars and text for tutorial purposes. The plot
shows the whole sample in grey and the chosen chromosome, not important in this guide, in red to blue gradient.
From this whole sample picture you can then use the arrangement of the clustered areas to estimate copy numbers
and allele-specific information. <br /><br />


 <img id="fig1taps" src="css/img/TAPS_call_help2.png" alt="Input arguments TAPS_call" title="Figure 1 TAPS"
 style="width:880px"><br /><br />

The <b>cn1</b>, <b>cn2</b> and <b>cn3</b> arguments are the positions of copy number 1, 2 and 3 on
the "Average log-ratio" axis. In this example cn2 is ~ -0.15. You need only provide any two of these values. <br />
The <b>loh</b> argument is the position of homozygous copy number 2 (LOH) on the "Allelic imbalance" axis.
In this example loh is ~ 0.75 

<br /><br /> -->

After running TAPS_call() you will have an additional set of plots in the sample folders which includes your
analysis of the plots.<br /><br />

In the results tab we will look at interpreting the data.

If you get any errors that you can't solve, please send us an 
<a onmouseover="popup('sebastian.dilorenzo(at)medsci.uu.se <br /> markus.mayrhofer(at)medsci.uu.se <br /> anders.isaksson(at)medsci.uu.se');">
email</a>
and we will get back to you as soon as possible!<br /><br />


<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>

