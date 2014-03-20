<!-- Flowchart -->
<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowchartTAPS.jpg" 
alt="Flowchart of processes and data for patchwork" title="Flowchart of data and processes for patchwork">
</div><br />

<center><h3>TAPS_plot()</h3></center>

TAPS_plot() is a function that visualizes allele-specific copy numbers from microarray
data of tumors. It can either be run in a directory containing multiple sample folders or
pointing to a specific sample, which is a folder containing no subfolders
and the files needed to execute.
<br /><br />

The execution time of TAPS_plot() varies. Samples that have been processed
by Nexus, where the segmentation step is already completed, may take anywhere from 5 minutes to half an hour. 
For other samples execution may take about an hour. <br /><br />

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

     TAPS_plot(directory=NULL,autoEstimate=FALSE,bin=400,cores=1,matched=FALSE)
     
Arguments:

directory: Default is NULL. If NULL getwd() is used. Specifying to a specific samples
          directory will run TAPS_plot on that directory.  Specifying a
          directory containing one or more subdirectories that are
          samples (and not any other subdirectories or TAPS_plot will
          error when trying to run them!) will iteratively run
          TAPS_plot on all samples.

autoEstimate: Default is FALSE. Development stage, do not use.

     bin: Default is 400.

   cores: Default is 1. Set amount of threads/cores to be used.

 matched: Default is FALSE. Set to TRUE if a matched normal has been
          used to remove homozygous snp allele frequencies.
</pre>

The "directory" argument can be set to a folder containing several sample folders, this will execute all sample folders. For example
setting directory to the Main_folder.<br /><br />

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

Here is an example execution, with output, of TAPS_plot() through TAPS_call() using the current working directory.
<br /><br />

<pre style="height:320px;overflow:scroll;">
> library(TAPS)

> TAPS_plot(cores=4)
[1] "Using working directory"
[1] "root: /home/User/Samples/"
[1] "1/4: Sample1 Loading"
[1] "2/4: Sample2 Loading"
[1] "3/4: Sample3 Loading"
[1] "4/4: Sample4 Loading"
[1] "3/4: Sample1 OK"
[1] "1/4: Sample3 OK"
[1] "4/4: Sample4 OK"
[1] "2/4: Sample2 OK"
>

> TAPS_estimates()
> 

> TAPS_click()
>

> TAPS_call(cores=4)
[1] "1/4: Sample1 Loading"
[1] "2/4: Sample2 Loading"
[1] "3/4: Sample3 Loading"
[1] "4/4: Sample4 Loading"
[1] "3/4: Sample1 OK"
[1] "1/4: Sample2 OK"
[1] "4/4: Sample4 OK"
[1] "2/4: Sample3 OK"
>
</pre>

This will have generated one "&lt;samplename&gt;_overview.jpeg" in the main directory (should
there be one), several plots and .Rdata files in the sample folder and a file named
simply "SampleData.csv" which will be covered shortly. <br /><br />

<hr class="alt1" />

<center><h3>User Evaluation</h3></center>

To infer the arguments for SampleData.csv that TAPS_call() uses you will need to look at one of
the chromosomal plots generated using TAPS_plot(). TAPS_estimate() and TAPS_click,(), discussed later, will help with the input
of the arguments. 
<br /><br />

The structure and relationships in the plot
can be interpreted to figure out the most probable locations of the allele-specific copy numbers.
For information to help you understand the structure of the plot and how the axis allelic imbalance
 and average log-ratio can be interpreted,
click <a href="AI_ALR.php" target="_blank" style="text-decoration:none;">here</a>.
<br /><br />

Now lets look at a general example of the structure and placement of clusters on the whole genome plot.<br /><br />
As average log-ratio increases, copy number increases. As Allelic imbalance increases, heterozygosity decreases,
that is to say that clusters with high allelic imbalance will be homozygous.<br /><br />
What do we expect a hypothetical plots arrangement of clusters to look like?
The average ploidy of the sample will be 0 on average log-ratio axis as it is normalized.
The sample may be highly rearranged but quite often this is a starting point for finding
copy number 2 or copy number 3. As there will be less reads covering copy number 1 in the sample
than higher copy numbers and copy number 1 cannot have different allele constitutions, by
its very nature of being one allele, copy number 1 will be represented by a single cluster far to the left
on average log-ratio axis when compared to the other clusters.<br /><br />
What if we do not have any copy number 1 in the sample? Then perhaps the far left
of the plot will be occupied by two clusters, indicating the homozygous, LOH, and heterozygous states
of copy number 2.
It then stands to reason that the next cluster we will encounter, again; moving from left to right
on average log-ratio axis, will be copy number 3.<br />
In the same way we would expect copy number 2 to follow
copy number 1 in the previous scenario. Note also that the lowest allelic imbalance of copy number 3
MUST be higher than the lowest allelic imbalance of copy number 2 as copy number 2 can have one of each allele
while copy number 3 must, at lowest allelic imbalance, have two of one allele and one of the other.<br />

The X chromosome is shown as "X"'s instead of gradient circles and may also be far left, do not let it
disturb your assessment of copy number levels as it can be misleading. <br /><br />

It is with reasoning such as this, looking at the plot and how the clusters are arranged and what cluster constitutions are
physically possible, that we can determine the allele-specific copy numbers in our sample.<br /><br />

Here is the topmost section of one such plot with some added bars and text for tutorial purposes. The plot
shows the whole sample in grey and the chosen chromosome, not important in this guide, in red to blue gradient.
From this whole sample picture you can then use the arrangement of the clustered areas to estimate copy numbers
and allele-specific information. <br /><br />


<img id="fig1taps" src="css/img/TAPS_call_help2.png" alt="Input arguments TAPS_call" title="Figure 1 TAPS"
 style="width:880px"><br /><br />

 The <b>cn1</b>, <b>cn2</b> and <b>cn3</b> arguments are the positions of copy number 1, 2 and 3 on
the "Average log-ratio" axis. In this example cn2 is ~ -0.15. You need only provide any two of these values later for TAPS_call(). <br />
The <b>loh</b> argument is the position of homozygous copy number 2 (LOH) on the "Allelic imbalance" axis.
In this example loh is ~ 0.75.
<br /><br />

We have developed two functions, TAPS_estimate() and TAPS_click(), that you can use if you do not want to fill SampleData.csv by hand.
This is especially useful if you are running multiple samples. We would recommend that you always use TAPS_click().
<br /><br />

<center><h4>TAPS_estimate()</h4></center>

TAPS_estimate() attempts to estimate the SampleData.csv arguments for all samples in the current working directory.
This information is hard to automate therefore it will not always do so correctly, but it will still simplify the process
when running TAPS_click() where you can either choose to accept TAPS_estimates() evaluation or redo it yourself.
<br /><br />

<center><h4>TAPS_click()</h4></center>
TAPS_click() does NOT work with Rstudio. R must be started from the commandline. <br /><br />

TAPS_click() launches a GUI for sample by sample interpretation of parameters. Much like TAPS_estimate(), the output of TAPS_click()
is in the form of inserting arguments for each sample into SampleData.csv. If TAPS_estimate() has been executed prior to
TAPS_click(), the first thing shown will be TAPS_estimates() interpretation, see plot 1 below.
At this stage you will be prompted by menu 1 and you will be given the choice to accept or redo the interpretation.
<br /><br />

<figure style="width:400px; margin:342px 0 0 0; float:left;">
<img src="css/img/TAPS_Menu1.png"><br /><br />
<i>Menu 1</i>, giving you choices to go to next sample (accept current estimates), redo the chromosome selection and interpretation of copy numbers of this sample, redo only
 the interpretation of copy numbers, skip the sample or close the program.
</figure>
<figure style="width:400px; margin:0 0 0 0; float:right;">
<img src="css/img/TAPS_Menu2.png"><br /><br />
<i>Menu 2</i>, choose to enter an interpretation of the copy numbers or look at a specific chromosome.
</figure>

<ul class="slideshow">
<li><img src="css/img/TAPS_Slide1.png" height="59%" width="59%" /><p>
<i>Plot 1</i>, TAPS_estimates() estimation of parameters for a sample. In this case the estimation is correct.
</p></li>
<li><img src="css/img/TAPS_Slide2.png" height="59%" width="59%" /><p>
  <i>Plot 2</i>, here we have chosen to perform copy number selection and have been prompted to identify via mouseclick the location of CN1.
</p></li>
<li><img src="css/img/TAPS_Slide3.png" height="59%" width="59%" /><p>
   <i>Plot 3</i>, TAPS_click() continues to prompt for input until completion. Here we have identified CN2.
</p></li>
<li><img src="css/img/TAPS_Slide4.png" height="59%" width="59%" /><p>
  <i>Plot 4</i>, CN3.
</p></li>
<li><img src="css/img/TAPS_Slide5.png" height="59%" width="59%" /><p>
  <i>Plot 5</i>, CN2LOH. The position of copy number 2 which has loss of heterozygosity. Located directly above copy number 2.
</p></li>
<li><img src="css/img/TAPS_Slide6.png" height="59%" width="59%" /><p>
  <i>Plot 6</i>, your new updated copy number interpretation. You will once again be prompted with menu 1 to move on to next sample etc.
</p></li>
</ul>


<hr class="alt1" />


<center><h3>TAPS_call()</h3></center>

The function TAPS_call() uses the relationship between log-ratio
and allelic imbalance to assign copy number and allele ratio for each segment.
<br /><br />

TAPS_plot() must have been executed before TAPS_call() can be. This is because the .Rdata files
generated during TAPS_plot() will be automatically read by TAPS_call() and you need to evaluate the plots
generated by TAPS_plot() to be able to give TAPS_call() meaningful input.  
<br /><br />

TAPS_call() has "directory" as its main argument, the parameters are read from the file
"SampleData.csv".
<br /><br />

An example call to TAPS_call(), where "mainsamplefolder" contains SampleData.csv: <br />

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

