TAPS is validated for output from the following arrays:<br />

<ul class="checks">
	<li>Affymetrix SNP6</li>
	<li>Affymetrix 250K/500K</li>
	<li>Affymetrix CytoScan HD</li>
</ul>

If you have some other data please contact us as it will likely work but as of writing is unvalidated.
<br /><br />

Your data can be prepared using Nexus (BioDiscovery, not free), ChAS (Affymetrix, free) or manually. 
Using Nexus has the benefit of using their segmentation which reduces running time of TAPS.
<br /><br />

If you have used <b>Nexus</b>; we recommend using the SNPRank segmentation algorithm as it is
CBS based rather than HMM, then use "File -> Utilities -> Export from .ivg to.txt" to create "probes.txt",
 "snps.txt" and "segments.txt". These should be in a sample folder, which is what TAPS takes as input.
<br /><br />

If you have used <b>ChAS</b>; load the file ending in "cyhd.cychp" into ChAS and create the text file ending in "cyhd.txt"
by choosing <b>Reports</b> menu and selecting <b>Export genotype results text file</b>. Then locate the samples folder,
which is what TAPS takes as input. It should now contain two files ending in "cyhd.txt" and "cyhd.cychp".
<br /><br />

If you want to create the required input files <b>manually</b>, store the Log-ratio as "probes.txt"
and allele frequency as "snps.txt" and place them in a folder, which will be your sample folder.
These two files should have this format to be compatible with TAPS: <br /><br />

<pre>
	$less probes.txt
		Chromosome      Start   End     Value  
		chr1    	15253   15278   0.2828165292739868      
		chr1    	48168   48193   0.2047460526227951      
		chr1    	60826   60851   -0.10802668333053589    
		chr1    	61722   61747   0.26339590549468994     
		chr1    	61795   61820   -0.1548314094543457     
		chr1    	61810   61835   0.11606350541114807
</pre>
<br />
<pre>
	$less snps.txt
		Chromosome      Start   End     Value
		chr1    	564621  564621  0.0     
		chr1    	721290  721290  0.0     
		chr1    	740857  740857  -5.551115123125783E-17  
		chr1    	752566  752566  0.29908037185668945     
		chr1    	761732  761732  1.1102230246251565E-16  
		chr1    	765269  765269  0.0     
</pre>
Both files are tab separated.

<br /><br /><div style="text-align:center;"><a href="#top" style="text-decoration:none;">Top</a></div>