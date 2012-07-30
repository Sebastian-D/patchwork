<!-- Flowchart -->
<!--<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowPWreq.png" 
alt="Flowchart of requirements for patchwork" title="Flowchart of requirements">
</div><br /> -->

To successfully run patchwork you will need to obtain/create these items:<br />

<ul class="checks">
	<!-- <li>A Human Genome fasta file</li> -->
	<li>A aligned and sorted BAM file of tumor content</li>
	<li>A BAI index of your BAM file</li>
	<li>A Pileup of your BAM file</li>
	<li>(optional) A matched normal sample to your tumor in BAM format</li>
	<li>(optional) A pileup of your normal sample BAM</li>
	<li>(optional) A BAI index of your normal sample BAM</li>
	<li>(optional) A standard Reference file. (Illumina/Solexa, SOLiD or your own)</li>
</ul>

<!--
The Human Genome fasta files can be found in many places. Here is a link for 1000genomes hosted
HG19 <a href="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
  target="_blank" style="text-decoration:none;">(human_g1k_v37)</a>. <br />
You probably also used such a file when aligning your tumor sequence.<br /><br /> -->

There are several (optional) files in the list, what we mean by this is that you do not need to have all these
files. There is really no point in using all of them.
<br /><br />
If you have a matched normal sample you should make a pileup of this and use those two arguments.
<br />
In some cases it may be better to use the reference file and a pileup.
<br />
You can also run patchwork.plot() with only a reference file.<br />
<br />

Here are two standard reference files for HG19. It is very important that you choose a reference
which matches the sequencing technology and reference genome used on your tumor sample.<br />

<a href="http://130.238.204.28/patchwork/references/datasolexa.RData" target="_blank" style="text-decoration:none;">
	Solexa/Illumina reference</a> <br />
<a href="http://130.238.204.28/patchwork/references/datasolid.RData" target="_blank" style="text-decoration:none;">
	SOLiD reference</a> <br /><br />

If you wish to create your own reference file, for example one that works for HG18, use patchworks.createreference().
For this you should use a pool of non-tumor BAM files, preferably no less than three, where the same sequencing
technology was applied as to your tumor sample. <br /><br />

<h4>patchwork.createreference()</h4> <br />

This function creates a reference using a pool of samples of your selection.
 Choosing the samples to use for reference creation you should consider these factors:

 <ul>
 	<li>Aligned using same reference genome (HG18/HG19)</li>
 	<li>They should be sequenced using the same technique</li>
 	<li>They should be from the same organism</li>
 	<li>They should be non-tumerous</li>
 </ul>

It is recommended that you use atleast 3 BAM files to create your reference.
<br /><br />

Initiate R and load the patchwork and patchworkData libraries: <br />

<pre>
	library(patchwork)
	library(patchworkData)
</pre>

Read the amazing documentation for patchwork.createreference():<br />

<pre>
	?patchwork.createreference
</pre>

The function takes a series of files with respective paths if they are not in your working directory, which is
where you started your R session by default. It also has a "output" parameter for customizing
output name. <br /><br />

Execute the function, pointing to your desired files: <br />

<pre>
	patchwork.createreference("file1","path/to/file2","../file3","~/heres/file4",output="REFOUT")
</pre>

This will generate REFOUT.Rdata, or whichever prefix you chose, in your working directory.
Use this file for the reference parameter of patchwork.plot(). <br /><br />

<h4>SAMtools file preparation</h4>

Many of the tasks that you will need to perform can be achieved using 
<a href="http://sourceforge.net/projects/samtools/files/" target="_blank" style="text-decoration:none;">
	SAMtools</a>. </br>
It is of special importance that you install SAMtools version 0.1.16 or older or the pileup generated will
have an incorrect format.<br />
An example of the correct format of your pileup file is included below. <br /><br />


To sort your BAM file: <br />


<pre>
	samtools sort &lt;tumorfile&gt;.bam &lt;sortedtumorfile&gt;.bam
</pre>

To create an index, BAI file, of either your normal sample or tumor sample BAM files: <br />

<pre>
	samtools index &lt;tumor/normalfile&gt;.bam
</pre>

The BAI file should always have the same name as your tumor file, for example "tumor.bam"
would have "tumor.bam.bai". If patchwork.plot() cannot find a correct BAI file this

<a onmouseover="popup('Traceback &#40;most recent call last&#41;&#58; <br /> File &quot;.python/Pysamloader.py&quot;, line 20, in &lt; module &gt; <p><dd> for read in infile.fetch&#40;a[1]&#41;&#58; </dd></p> File &quot;csamtools.pyx&quot;, line 745, in csamtools.Samfile.fetch &#40;pysam/csamtools.c&#58;8332&#41;<br />File &quot;csamtools.pyx&quot;, line 684, in csamtools.Samfile._parseRegion &#40;pysam/csamtools.c&#58;7615&#41; <br />ValueError: invalid reference &lsquo;1&lsquo;<br />Error in read.table(file = file, header = header, sep = sep, quote = quote,  &#58; <p><dd> no lines available in input </p></dd> Calls&#58; patchwork.plot -&gt; patchwork.readChroms -&gt; read.csv -&gt; read.table');">
error</a> is displayed.<br /><br />


To pileup either your normal sample or tumor sample BAM files: <br />

<pre>
	samtools pileup -vcf &lt;humangenome&gt;.fasta &lt;tumor/normalfile&gt;.bam &gt; pileup
</pre>

Your pileup file should have this format: <br />

<pre>
	$less pileup
		chr1    10179   c       W       0       0       60      5       ,,tna   ##6!3
		chr1    10180   t       W       6       6       60      5       ,,,aa   ##369
		chr1    10377   a       R       0       1       60      1       g       @
		chr1    11391   t       A       0       3       60      1       a       B
		chr1    18592   C       Y       0       3       60      1       t       B
		chr1    23359   c       S       0       2       60      1       g       A
		chr1    24067   A       R       0       2       60      1       G       A
		chr1    30315   G       C       3       24      60      2       cC      =;
		chr1    92200   a       W       0       3       60      1       t       B
		chr1    96592   T       C       0       3       60      1       C       B
		chr1    100140  a       M       0       1       60      1       C       @
		chr1    104697  g       K       0       2       60      1       T       A
		chr1    127285  A       R       0       3       60      1       g       B
</pre>

You should now have all the components needed for execution!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>

