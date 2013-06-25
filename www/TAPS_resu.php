<h4>Plots of TAPS_plot()</h4>

After running TAPS_plot() you should have chromosomal plots and an overview plot in your sample folder. 
<br /><br />

For information regarding allelic imbalance and average log-ratio click <a href="AI_ALR.php" target="_blank" style="text-decoration:none;">here</a>.
<br /><br />

<h4>Overview plot</h4>

<img id="fig2taps" src="css/img/TAPS_overview2.jpg" alt="Results example TAPS_plot() overview"
 title="Figure 2 TAPS"
 style="width:880px;"><br /><br />

<h5>Top half of the plot</h5>

All the possible chromosomes plotted for allelic imbalance versus average log-ratio.
The colored segments represent a specific chromosome while the grey is the rest of the sample.
If you look closely there are three grey bars in each plot.
The middle bar is located at 0 average log-ratio, the ploidy of the sample.

<h5>Bottom half of the plot</h5>

Log-ratio, cytoband and allele frequency information for each chromosome. The colored segments in log-ratio correspond
to the chromosome specific colored segments seen in the top part of the plot.
<br /><br />

<hr class="alt1" />
<h4>Chromosomal plots</h4>

<img id="fig3taps" src="css/img/TAPS_chromosomal2.jpg" alt="Results example TAPS_plot() chromosomal"
 title="Figure 3 TAPS"
 style="width:880px;"><br /><br />

<h5>Top half of the plot</h5>

The clusters in the plot display regions of certain allelic constitution and copy number.
The copy number increases along the average log-ratio axis while paternal/maternal allele ratio
becomes less balanced along the allelic imbalance axis.<br /><br />

The chromosome in question is colored against a background of the rest of the sample in grey.
A colored circles gradient correlates with its segments position, seen in log-ratio, on the
chromosome.<br />
The circles are semi-transparent so a darker hue, both for colored and grey,
indicate a greater amount of genomic content in that region.<br /><br />

The interpretation of this part of the plot especially in allele-specific copy numbers is
the input for TAPS_call() and has been stated in TAPS execution tab.<br /><br />

Quoted here for your convenience:<br />
<adress><p>
	"Lets take a look at the structure and placement of clusters on the whole genome plot.<br /><br />
	As average log-ratio increases, copy numbers increases. As allelic imbalance increases, heterozygosity decreases,
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
	physically possible, that we can determine the allele-specific copynumbers in our sample."<br /><br />
</p></adress>


<h5>Bottom half of the plot</h5>
Log-ratio, cytoband and allele frequency information for the selected chromosome.
The colored segments in log-ratio correspond directly to the colored segments seen in the top part of the plot.
<br /><br />
Higher copy numbers will have a higher log-ratio.
<br /><br />
Allele frequencies detected will vary depending on the allelic constitution of a segment of a chromosome.
Seen here in the very frist part of the chromosome, between 0-30mb approximately, where the sample is copy number
2 with loss of heterozygosity. This means there is no variation in the allele constitution for that segment.
Looking on at the rest of the chromosome in question it is copy number 3 with 2 of one allele and 1 of the other.
This is shown in the banded structure, given by presence of more than one option, of allele frequency as well.
<br /><br />

<hr class="alt1" />

<h4>Plots of TAPS_call()</h4>

TAPS_call() generates you interpretation, as specified in execution tab, of the allele-specific copy numbers in the sample.
It generates chromosomal plots similar to the ones of TAPS_plot() but additionally specifies total and minor
copy numbers of segments. It also generates a quick "check" overview of the sample with all allele-specific copy numbers
specified.
<br /><br />

<h4>Check</h4>

<img id="fig4taps" src="css/img/TAPS_check2.png" alt="Results example TAPS_call() check"
 title="Figure 4 TAPS"
 style="width:880px;"><br /><br />

This plot, arbitrarily named "check", is not chromosome specific and shows
the complete tumor sample with assigned tags for copy number and allele content,
much like the <a href="css/img/TAPS_AI_ALR.jpeg" target="_blank" style="text-decoration:none;">explanatory picture</a>
seen previously but with text describing the copy numbers and allele ratios of the segments
rather than purple and green bars. 
<br /><br />

<b>1m0</b> is copy number 1, homozygous. <br />
<b>2m1</b> is copy number 2, heterozygous. <br />
<b>2m0</b> is copy number 2, homozygous. <br />
<b>3m1</b> is copy number 3, heterozygous. <br />
<b>3m0</b> is copy number 3, homozygous. <br />
And so on. <br /><br />

If the texts showing the approximate positions of the allele-specific copy numbers seem missplaced you may want
to try to improve, or tweak, your input values in SampleData.txt, the parameters used by TAPS_call().
<br /><br />

<hr class="alt1" />

<h4>Chromosomal plots with total and minor copynumbers</h4>

The "by chromosome" plots generated by TAPS_call() are quite similar to the chromosomal plots generated by
TAPS_plot() with the addition that they now include your allele-specific copy number analysis to display
total and minor copy numbers.

<img id="fig5taps" src="css/img/TAPS_chromosomalCN.jpg" alt="Results example TAPS_call() chromosomal total and minor copy numbers"
 title="Figure 5 TAPS"
 style="width:880px;"><br /><br />

<h5>Left half of the plot</h5>

The whole sample with the current chromosomes, in this example chromosome 7, segments colored from blue to red.
The X chromosome is shown as X's. The size of the circles vary with length of segment, larger for larger segments.
Allele-specific copy numbers of the different clusters, 1m0, 2m1 etc, are shown for quick interpretation.

<h5>Right half of the plot</h5>

Farthest up we see the total copy number, colored, and minor copy number, grey, of a segment.
If a segments total copy number is 3 and 
minor copy number is 1 then that segment of the chromosome has 1 paternal and 2 maternal
OR 2 paternal and 1 maternal copies.<br />
Another example; if a segments total copy number is 2 and minor copy number is 0 then
that segment of the chromosome has 2 paternal OR 2 maternal copies.
<br /><br />

Below that we see the same three plots of log-ratio, cytoband information and allele frequency as seen above
in chromosomal TAPS_plot() plots. Notice that all four of the plots on the right side have the position on the
chromsome in question as their x axis.

<br /><br />For more background information please contact us!

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>



















