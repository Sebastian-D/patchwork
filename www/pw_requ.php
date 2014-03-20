<!-- Flowchart -->
<!--<div style="text-align:center;">
<img id="flowchartcg" src="css/img/flowPWreq.png" 
alt="Flowchart of requirements for patchwork" title="Flowchart of requirements">
</div><br /> -->

Here is information on the files you need and how to get them.<br /><br />
To successfully run patchwork you will need to obtain/create these items:<br />

<ul class="checks">
	<!-- <li>A Human Genome fasta file</li> -->
	<!-- <li>An aligned and sorted BAM file of tumor content</li>
	<li>A BAI index of your BAM file</li> -->
	<li><b>An aligned and sorted tumor BAM file.</b> <br />A BAI index of your BAM file. <br />A pileup of your BAM file. <br />A VCF file of your tumor pileup (if applicable).</li>
	<!-- <li>(optional) A VCF file of your tumor pileup</li> -->
	<li><b>(optional) A matched normal sample to your tumor in BAM format.</b><br /> A BAI index of your normal sample BAM.<br /> A pileup of your normal sample BAM. <br />A VCF file of your normal pileup (if applicable). </li>
	<!-- <li>(optional) A pileup of your normal sample BAM</li>
	<li>(optional) A VCF file of your normal pileup</li>
	<li>(optional) A BAI index of your normal sample BAM</li> -->
	<li><b>(optional) A standard Reference file. (Illumina/Solexa, SOLiD or your own)</b></li>
</ul>

<!--
The Human Genome fasta files can be found in many places. Here is a link for 1000genomes hosted
HG19 <a href="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
  target="_blank" style="text-decoration:none;">(human_g1k_v37)</a>. <br />
You probably also used such a file when aligning your tumor sequence.<br /><br /> -->

The (optional) files are interchangeable. If you have a normal sample you do not need to supply a reference file.
<br /><br />
If you have a matched normal sample you should create a pileup and BAI file.
The BAI file is required for the file to be read by patchwork.
<br />
You can also use a reference file, see below, if you do not have a matched normal.<br />
<br />

The optional VCF files need to be supplied if you are using a later version of samtools and have opted
for using samtools mpileup to create the pileup file. The vcfs contain consensus and variant information which
the mpileup generated pileup lacks. We recommend using the pileup from older version of samtools rather than mpileup + vcf of new version.
More information on generating these files is further down.<br /><br />

Here are two standard reference files created for patchwork for samples where alignment used the UCSC HG19 reference.
It is very important that you choose a reference which matches the sequencing technology and reference genome used for
 your tumor sample.<br />

<a href="http://array.medsci.uu.se/patchwork/datasolexa.RData" target="_blank" style="text-decoration:none;">
	Solexa/Illumina reference</a> <br />
<a href="http://array.medsci.uu.se/patchwork/datasolid.RData" target="_blank" style="text-decoration:none;">
	SOLiD reference</a> <br /><br />

If you wish to create your own reference file, for example one that works for HG18, use patchworks.createreference().
<!--For this you should use a pool of non-tumor BAM files, preferably no less than three, where the same sequencing
technology was applied as to your tumor sample. --> <br /><br />

<center><h3>patchwork.createreference()</h3></center> <br />

This function creates a reference using a pool of samples of your selection.
 Choosing the samples to use for reference creation you should consider these factors:

 <ul>
 	<li>Aligned using same reference genome (HG18/HG19)</li>
 	<li>Sequenced using the same technique</li>
 	<li>From the same organism</li>
 	<li>Non-tumorous</li>
 	<li>Same gender</li>
 </ul>

It is recommended that you use at least 3 BAM files to create your reference.
Note also that if your tumor sample is a different gender than your reference samples you may not get optimal results.
There is a paramter in the second part of analysis, patchwork.copynumbers(), to correct for this. If you are specifically interested
in the sex chromosomes we would definitely recommend using the same gender in reference as in tumor sample.
<br /><br />

To execute patchwork.createreference() start R and load the patchwork and patchworkData libraries: <br />

<pre>
	library(patchwork)
	library(patchworkData)
</pre>

Read the documentation for patchwork.createreference():<br />

<pre>
	?patchwork.createreference
</pre>

The function takes a series of BAM files as input and outputs a single reference file. 
It also has a "output" parameter for customizing output name. <br /><br />

Execute the function, pointing to your desired files: <br />

<pre>
	patchwork.createreference("file1","path/to/file2","../file3","~/here_is/file4",output="REFOUT")
</pre>

This will generate REFOUT.Rdata, or whichever prefix you chose, in your working directory.
Use this file for the reference parameter of patchwork.plot(). <br /><br />

<center><h3>SAMtools file preparation</h3></center>

Many of the tasks that you will need to perform can be achieved using 
<a href="http://sourceforge.net/projects/samtools/files/" target="_blank" style="text-decoration:none;">
	SAMtools</a>. </br>
If you wish to save space and computation time you should install SAMtools version 0.1.16 or older, otherwise the pileup generated will
have an incorrect format. Instructions for using the newest version of SAMtools or older version is below.<br /><br />

To sort your BAM file: <br />

<pre>
	samtools sort &lt;tumorfile&gt;.bam &lt;sortedtumorfile&gt;.bam
</pre>

To create an index, BAI file, of either your normal sample or tumor sample BAM files: <br />

<pre>
	samtools index &lt;tumor_or_normalfile&gt;.bam
</pre>

The BAI file should always have the same name as your tumor file, for example "tumor.bam"
would have "tumor.bam.bai". If patchwork.plot() cannot find a correct BAI file this

<a onmouseover="popup('Traceback &#40;most recent call last&#41;&#58; <br /> File &quot;.python/Pysamloader.py&quot;, line 20, in &lt; module &gt; <p><dd> for read in infile.fetch&#40;a[1]&#41;&#58; </dd></p> File &quot;csamtools.pyx&quot;, line 745, in csamtools.Samfile.fetch &#40;pysam/csamtools.c&#58;8332&#41;<br />File &quot;csamtools.pyx&quot;, line 684, in csamtools.Samfile._parseRegion &#40;pysam/csamtools.c&#58;7615&#41; <br />ValueError: invalid reference &lsquo;1&lsquo;<br />Error in read.table(file = file, header = header, sep = sep, quote = quote,  &#58; <p><dd> no lines available in input </p></dd> Calls&#58; patchwork.plot -&gt; patchwork.readChroms -&gt; read.csv -&gt; read.table');">
error</a> is displayed.<br /><br />

<center><h4>Instructions for SAMtools &lt;= 0.1.16 </h4></center>

To pileup either your normal sample or tumor sample BAM files: <br />

<pre>
	samtools pileup -vcf &lt;humangenome&gt;.fasta &lt;tumor_or_normalfile&gt;.bam &gt; pileup
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

<center><h4>Instructions for SAMtools &gt; 0.1.16 </h4></center>

To mpileup either your normal sample or tumor sample BAM files: <br />

<pre>
	samtools mpileup -f &lt;humangenome&gt;.fasta &lt;tumor_or_normal&gt;.bam &gt; mpileup
</pre>

For consensus calling: <br />

<pre>
	samtools mpileup -uf &lt;humangenome&gt;.fasta &lt;tumor_or_normal&gt;.bam |
		bcftools view -bvcg - &gt; &lt;unfiltered_output&gt;.bcf

	bcftools view &lt;unfiltered_output&gt;.bcf | vcfutils.pl varFilter -D100 > &lt;output&gt;.vcf
</pre>

The -D100 option filters out SNPs with a read depth higher than 100
(They might be from repeat regions in your sample), you should change this parameter based on the average coverage of
your sample.<br />
For average coverage 30x =~ -D100, average coverage 120x =~ -D500 etc. <br /><br />

Remember to use both of these files when running patchwork.plot() with the vcf file in the Tumor.vcf/Normal.vcf parameter
 and the mpileup in Tumor.pileup/Normal.pileup parameter. <br /><br />

You should now have all the components needed for execution!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>

