<html><head>

<script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.6.4/jquery.min.js"></script>
<script type="text/javascript" src="js/prettify.js"></script>
                                <!-- PRETTIFY -->
<script type="text/javascript" src="js/kickstart.js"></script>
<link rel="stylesheet" type="text/css" href="css/kickstart.css" media="all" />
<link rel="stylesheet" type="text/css" href="style.css" media="all" />

</head>
<body>

<div id="wrap" class="clearfix">
	<div class="col_12">

		<h4>Allelic Imbalance</h4>

		Allele-specific information is obtained from read counts at heterozygous SNPs.

		<!--
		Allelic imbalance of a segment is then calculated as 


		<table class="fraction" cellpadding="0" cellspacing="0">
			<tr>
				<td style="text-align:center;" nowrap="nowrap"><i>&Sigma;(Highest allele counts)
				&minus;&Sigma;(Lowest allele counts)</i>
				</td>
			</tr>
			<tr>
				<td style="text-align:center;" class="upper_line"><i>&Sigma;(Highest allele counts)</i>
				</td>
			</tr>
		</table>


		What this basically does is separate chromosomal regions by their allele constitution. -->
		Allelic-imbalance separates chromosomal regions by their allele constitution.
		A diploid cell, which is the norm, will have one maternal and one paternal copy of a chromosome.
		For each copy number there can be losses of segments to either of the
		chromosomes as well as copies/gains. It is this allele-specific information that
		allelic imbalance captures. It does not, however, discern the origin of a segment,
		that is to say that two paternal copies have the same allelic imbalance ratio as two
		 maternal copies. Below is an attempt to explain the different possible allelic combinations
		 for a few copy numbers,
		 it is VERY important that you understand that it is the ratio between green and purple that
		 is important, not if purple/green is maternal/paternal. <br />
		 From an allelic imbalance score point of view
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;"> is the same as 
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">. <br /><br />

		<table border="1">

		<tr style="height:66px">
			<th style="border-bottom:1px solid black;">Copy Number</th>
			<th colspan=3 style="border-bottom:1px solid black;">Possible chromosomal segment combinations</th>
		</tr>

		<tr>
			<td>1</td>
			<td><img src="css/img/purple_chromosome.png" alt="purple segment"
			title="Purple segment" style="display:inline-block;"></td>
		</tr>

		<tr>
			<td>2</td>
			<td><img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td><img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
		</tr>

		<tr>
			<td>3</td>
			<td>
			<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td><img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
		</tr>

		<tr>
			<td>4</td>
			<td><img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td><img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td>
			<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
		</tr>

		<tr>
			<td>5</td>
			<td><img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td>
			<img src="css/img/purple_chromosome.png" alt="purple segment"
		title="Purple segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
			<td>
			<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;">
		<img src="css/img/green_chromosome.png" alt="green segment"
		title="Green segment" style="display:inline-block;"></td>
		</tr>

		</table><br />

		Here is a an example of a plot generated using patchwork where I have added the possible
		allele combinations next to the clusters. <br />

		<img src="css/img/Allelic_plot.png" alt="Allele-specific tutorial" title="Allele-specific graphic guide"
		 style="width:880px;">

		<br /><h4>Coverage</h4>

		Coverage is measured as total DNA content and normalized around 1. You will get a 
		clearer picture at higher coverage data, 120x etc, however patchworkCG has high quality results
		with low coverage aswell.

	</div>
	<div class="clear"></div>
</div>

<div id="footer">
	<div id="left">
		Sebastian DiLorenzo 	-	sebastian.dilorenzo(at)medsci.uu.se <br />

	   	Markus Mayrhofer		-	markus.mayrhofer(at)medsci.uu.se <br />

	   	Anders Isaksson		-	anders.isaksson(at)medsci.uu.se <br />

		   	<br />
		   	<br />
		   	Last modified: 2012-04-16

	</div>
	<div style="text-align:right;" id="right">
		<a href="http://r-forge.r-project.org/" target="_blank" style="text-decoration:none;">
			R-forge</a> <br />
	    <a href="http://www.completegenomics.com" target="_blank" style="text-decoration:none;">
	    	CompleteGenomics</a> <br />
	    <a href="http://www.scilifelab.se/index.php?content=home" target="_blank" style="text-decoration:none;">
	    	SciLifeLab</a><br />
	    <a href="http://www.akademiska.se/" target="_blank" style="text-decoration:none;">
	    	Akademiska Sjukhuset</a> <br />
	    <a href="http://www.medsci.uu.se" target="_blank" style="text-decoration:none;">
	    	Uppsala University Medical Sciences</a> <br />
	</div>
	<div class="clear"></div>
</div>
		
</body>
</html>











