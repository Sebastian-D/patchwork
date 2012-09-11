<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">
<html><head>
<title>Patchwork</title>
<meta charset="UTF-8">

<link type="text/plain" rel="author" href="http://patchwork.r-forge.r-project.org/humans.txt" />

<meta name="description" content="Tutorial and information regarding Patchwork and PatchworkCG, an allele-specific visualization tool developed for cancer genomes." />
<meta name="author" content="Sebastian DiLorenzo" />
<meta name="keywords" content="Cancer,patchwork,patchworkCG,CompleteGenomics,R,visualization,Sequencing,Rforge,NGS,WGS,Allelic,Imbalance,Coverage,Allele-specific,Copy number" />

<script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.6.4/jquery.min.js"></script>
<noscript>
Your browser does not support JavaScript! Enable it for this page to work as intended.
</noscript>
<!-- <script type='text/javascript' src='js/jquery.ba-hashchange.min.js'></script>
<script type='text/javascript' src='js/dynamicpage.js'></script> -->
<!--[if lt IE 9]><script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script><![endif]-->
<script type="text/javascript" src="js/prettify.js"></script>
                                <!-- PRETTIFY -->
<script type="text/javascript" src="js/kickstart.js"></script>
<link rel="stylesheet" type="text/css" href="css/kickstart.css" media="all" />
                     	<!-- CUSTOM STYLES -->
<link rel="stylesheet" type="text/css" href="style.css" media="all" />    


<script type="text/javascript" src="js/sidebar.js"></script>

<script type="text/javascript" src="js/popup.js"></script>

<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-31529932-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

</script>

</head>

<body>

<a id="top-of-page"></a>
<a name="top"></a>

<!-- collapsible sidebar -->
<!-- <div id="colleft">
<div id="panel" style="border:1px solid #CCC;">
	<a href="http://www.medsci.uu.se" accesskey="1" tabindex="2">
			<img  height="118" width="118" 
			src="http://www.medsci.uu.se/digitalAssets/22/22975_UU_logga_transp.png" 
			alt="Medical Sciences" title="Medical Sciences Uppsala University" >
	</a>
	<a href="http://www.completegenomics.com/">
	   <img src="http://cgatools.sourceforge.net/cg-logo.png" height="36px" width="116px" title="Complete Genomics"
	   alt="Complete Genomics homepage">
	 </a>

	 <a href="http://r-forge.r-project.org/">
		<img src="http://r-forge.r-project.org/themes/rforge/imagesrf/logo.png" style="padding-top:10px;" 
		height="36" width="116" alt="R-Forge Logo" title="R-Forge">
	</a>

	<div style="text-align:center;">
	<a href="http://www.akademiska.se">
		<img src="http://www.akademiska.se/upload/grafisk_profil/jpg/34mm/aka_34mm_72dpi.jpg"
		title="Akademiska Sjukhuset" alt="Akademiska Hospital" style="padding-top:10px;">
	</a>
	</div>

	<a href="http://www.scilifelab.se">
		<img src="http://www.scilifelab.se/images/logo_header.png"
		title="SciLifeLab" alt="SciLifeLab Homepage"
		height="36" width="116" style="padding-top:10px;">
	</a>

<div id="hidePanel" style="padding-top:10px;"><a href="#" style="text-decoration:none; padding-left:40px;">
<span style="color:#999; font-size:24px;">
 &laquo;</span></a></div>
</div></div>
<div id="showPanel"><span style="color:#999; border:1px solid #CCC;">&raquo;</span></div>
-->

<div id="wrap" class="clearfix">
<!-- ======= END HEADER ====== -->

<!-- <div id="container"> -->
	<div id="left">
 		<h1 style="padding-left:20px; color:#C7C7C7;"><a href="http://patchwork.r-forge.r-project.org/" style="text-decoration:none; color:#C7C7C7;">Patchwork</a></h1>
 	</div>

<!-- 	<div id="right"> -->
<!-- 		<a href="http://www.medsci.uu.se" accesskey="1" tabindex="2">
			<img class="align-right" height="125" width="125" 
			src="http://www.medsci.uu.se/digitalAssets/22/22975_UU_logga_transp.png" 
			alt="Medical Sciences" title="Medical Sciences Uppsala University" >
		</a> -->
<!-- </div> -->
<!-- </div> -->

<!--	
	<a href="http://www.completegenomics.com/">
	   <img class="align-right" id="img_logo" src="http://cgatools.sourceforge.net/cg-logo.png" height="53px" width="210px" style="padding-top:20px" alt="Complete Genomics">
	 </a>

	<a href="http://r-forge.r-project.org/">
		<img class="align-right" src="http://r-forge.r-project.org/themes/rforge/imagesrf/logo.png" style="padding-top:20px" alt="R-Forge Logo" />
	</a>
-->
 

	<div class="col_12">
	<!-- Tabs Left -->
		<ul class="tabs left">
			<li><a href="#tabr0">Introduction</a></li>
			<li><a href="#tabr1">PatchworkCG</a></li>
			<li><a href="#tabr2">Patchwork</a></li>
			<li style="margin:0 0 0 52%;"><a href="#tabr3">Changelog</a></li>
		</ul>

<!-- HOME/INTRODUCTION -->
		<div id="tabr0" class="tab-content">

		<?php include("intro.php"); ?>

		</div>

<!-- PATCHWORKCG -->
		<div id="tabr1" class="tab-content" style="padding:0; margin:0;">
				<ul class="tabs right">
					<li><a href="#tabr20"><span class="icon" data-icon="g"></span>Installation</a></li>
					<li><a href="#tabr21"><span class="icon" data-icon=","></span>Requirements</a></li>
					<li><a href="#tabr22"><span class="icon" data-icon="V"></span>Execution</a></li>
					<li style="margin:0 133px 0 0;"><a href="#tabr23"><span class="icon" data-icon="j"></span>Results</a></li>
				<li><a href="#tabr24"><span class="icon" data-icon="i"></span>Demo Data</a></li>
				</ul>

			<div id="tabr20" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("cg_inst.php"); ?>
			</div>
			<div id="tabr21" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("cg_requ.php"); ?>
			</div>
			<div id="tabr22" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("cg_exec.php"); ?>
			</div>
			<div id="tabr23" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("cg_resu.php"); ?>
			</div>
			<div id="tabr24" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("cg_demo.php"); ?>
			</div>
		</div>

<!-- PATCHWORK -->
		<div id="tabr2" class="tab-content" style="padding:0; margin:0;">
			<ul class="tabs right">
				<li><a href="#tabr10"><span class="icon" data-icon="g"></span>Installation</a></li>
				<li><a href="#tabr11"><span class="icon" data-icon=","></span>Requirements</a></li>
				<li><a href="#tabr12"><span class="icon" data-icon="V"></span>Execution</a></li>
				<li style="margin:0 133px 0 0;"><a href="#tabr13"><span class="icon" data-icon="j"></span>Results</a></li>
				<li style="visibility:hidden;"><a href="#tabr14"><span class="icon" data-icon="i"></span>Demo Data</a></li>
			</ul>

			<div id="tabr10" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("pw_inst.php"); ?>
			</div>
			<div id="tabr11" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("pw_requ.php"); ?>
			</div>
			<div id="tabr12" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("pw_exec.php"); ?>
			</div>
			<div id="tabr13" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("pw_resu.php"); ?>
			</div>
			<!-- <div id="tabr14" class="tab-content" style="border:1px solid transparent; border-left:1px solid transparent;">
				<?php include("pw_demo.php"); ?>
			</div> -->
		</div>

		<!-- CHANGELOG -->
		<div id="tabr3" class="tab-content">
			<?php include("changelog.php"); ?>
		</div>

	</div>
	<div class="clear"></div>
</div><!-- END WRAP -->

<!-- ============ START FOOTER ================ -->
  
	<div id="footer">
		<div id="left">
			Sebastian DiLorenzo - sebastian.dilorenzo(at)medsci.uu.se <br />

		   	Markus Mayrhofer - markus.mayrhofer(at)medsci.uu.se <br />

		   	Anders Isaksson	- anders.isaksson(at)medsci.uu.se <br />

		   	<br />
		   	<br />
		   	
		   	<!-- Start of StatCounter Code for Default Guide -->
			<script type="text/javascript">
				var sc_project=7876707; 
				var sc_invisible=0; 
				var sc_security="39e2d4af"; 
			</script>
			<script type="text/javascript"
				src="http://www.statcounter.com/counter/counter.js">
			</script>
			<noscript>
				<div class="statcounter">
					<a title="drupal
						statistics module" href="http://statcounter.com/drupal/"
						target="_blank"><img class="statcounter"
						src="http://c.statcounter.com/7876707/0/39e2d4af/0/"
						alt="drupal statistics module">
					</a>
				</div>
			</noscript>
			<!-- End of StatCounter Code for Default Guide -->
			<div style="display:inline-block; width:2px; height:1px;"></div>
			<!-- LAST MODIFIED -->
			Last modified: 2012-09-11


		</div>
		<div style="text-align:right;" id="right">
			<a href="http://r-forge.r-project.org/" target="_blank" style="text-decoration:none;">
				R-forge</a> <br />
		    <a href="http://www.completegenomics.com" target="_blank" style="text-decoration:none;">
		    	CompleteGenomics</a> <br />
		    <a href="http://www.scilifelab.se/index.php?content=home" target="_blank" style="text-decoration:none;">
		    	SciLifeLab</a> <br />
		    <a href="http://www.akademiska.se/" target="_blank" style="text-decoration:none;">
		    	Akademiska Sjukhuset</a> <br />
		    <a href="http://www.medsci.uu.se" target="_blank" style="text-decoration:none;">
		    	Uppsala University Medical Sciences</a> <br />
		</div>
		<div class="clear"></div>
		<!-- </div> -->
	</div>



</body>
</html>