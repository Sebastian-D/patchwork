<h4>Plots of patchwork.plot()</h4>

After running patchwork.plot() you should have 24 plots in your working directory similar
to the example plot below, one for each chromosome.<br /><br />

For information regarding allelic imbalance and coverage click <a href="AI_Cov.php" target="_blank" style="text-decoration:none;">here</a>.

<h5>Top part of the plot</h5>

The clusters in the plot display regions of certain allelic constitution and copy number.
The copy number increases along the Coverage axis while paternal/maternal allele ratio
becomes less balanced along the Allelic Imbalance axis.<br /><br />

The chromosome in question is colored against a background of the complete genome in grey.
A colored circles gradient and size correlate with its segments position and size on the
 chromosome. Colors of segments match between plots.
 The circles are semi-transparent so a darker hue, both for colored and grey,
 indicate a greater amount of genomic content in that region.<br /><br />

We know that each cluster has a certain copy number and allele content and we know that
a human is usually diploid (copy number 2, heterozygous). Finally we know that the average
copy number of the genome in question is at position 1 on the Coverage axis.
We can use this information to determine the clusters probable copy number
and allele content.<br /><br />

The far left tiny cluster, on Coverage axis, is the deletions, copy number 1.
Moving to the right the next two clusters, lower cluster is quite small,
are the allelic states of copy number 2. Following this reasoning the next set of clusters
must be copy number 3, then 4 etc. The average copy number for the sample is just above 3.

This arrangement of the genome is easier to see from one of the plots generated by
patchwork.copynumbers() or in our section above on allelic imbalance, however before using
 that function you will need to be able to determine the constitution of the tumor
  genome as it is used as input arguments for patchwork.copynumbers()!

<img id="fig2pw" src="css/img/tumor.bam_karyotype.chr18.png" alt="Results example patchwork.plot"
 title="Figure 2 patchwork"
 style="width:880px;"><br /><br />


 <h5>Middle part of the plot</h5>

The chromosome in questions coverage plotted against the chromosomes coordinates.

<h5>Bottom part of the plot</h5>

The chromosome in questions allelic imbalance plotted against the
 chromosomes coordinates.

<br /><br /><h4>Plots of patchwork.copynumbers()</h4>

This plot, arbitrarily named "check", is not chromosome specific and shows
the complete tumor genome with assigned tags for copy number and allele content,
much like the <a href="css/img/Allelic_plot.png" target="_blank" style="text-decoration:none;">explanatory picture</a>
seen previously but with text describing the copy numbers and allele ratios of the segments
rather than purple and green bars.

<img id="fig3pw" src="css/img/copynumbers__karyotype_check.png" alt="Results example patchwork check"
 title="Figure 3 patchwork"
 style="width:880px"><br /><br />

<b>1m0</b> is copy number 1, homozygous. <br />
<b>2m1</b> is copy number 2, heterozygous. <br />
<b>2m0</b> is copy number 2, homozygous. <br />
<b>3m1</b> is copy number 3, heterozygous. <br />
<b>3m0</b> is copy number 3, homozygous. <br />
And so on. <br /><br />

The other chromosomal plots generated from patchwork.copynumbers() are the same as
the ones generated by patchwork.plot() but with one extra part of the plot. <br />

<!-- change picture -->
<img id="fig4pw" src="css/img/copynumbers__karyotypeCN.chr18.png" alt="Results example patchwork.copynumbers"
 title="Figure 4 patchwork"
 style="width:880px"><br /><br />

In the part of the plot seen second from the top we see the total copy number, colored,
and minor copy number, grey, of a segment. If a segments total copy number is 3 and 
minor copy number is 1 then that segment of the chromosome has 1 paternal and 3 maternal
OR 3 paternal and 1 maternal copies.<br />
Another example; if a segments total copy number is 2 and minor copy number is 0 then
that segment of the chromosome has 2 paternal OR 2 maternal copies.

<br /><br />For more background information please read our submitted article!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>


































