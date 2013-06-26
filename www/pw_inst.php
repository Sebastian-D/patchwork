Patchwork works on Unix (Linux, Ubuntu, MacOSX, etc) based systems. You will probably get an error
when trying to install it on Windows, and you definitely can't run it.
<br /><br />

Begin by starting R. It is recommended that you use the latest version.

<br /><br />
To install DNACopy, which is a program patchwork needs to run, issue the following commands within R.

<pre>
    source("http://bioconductor.org/biocLite.R")
    biocLite("DNAcopy")
</pre>

After DNAcopy is installed you can install patchwork and patchworkData.

<pre>
    install.packages("patchworkData", repos="http://R-Forge.R-project.org")
    install.packages("patchwork", repos="http://R-Forge.R-project.org")
</pre>

If for some reason that does not work add the 'type="source"' to it, as so:

<pre>
    install.packages("patchworkData", repos="http://R-Forge.R-project.org",type="source")
    install.packages("patchwork", repos="http://R-Forge.R-project.org",type="source")
</pre>

If something still does not work with installation, here are links to the location of the source files: <br />

<a href="http://r-forge.r-project.org/R/?group_id=1250" target="_blank" style="text-decoration:none;">
	Patchwork</a><br />
<a href="http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html" target="_blank" style="text-decoration:none;">
	DNACopy</a>

<br />

<!--
You only need to install patchwork once. But remember to keep it up to date!

<pre>
    update.packages(repos=\"http://R-Forge.R-project.org\")
</pre>
-->

<h4>Python and Pysam</h4>

<!-- Patchwork uses a python mod called pysam for "reading" the chromosomes. The package already checks if you have
it installed and otherwise installs it for you the first time you run patchwork however if you lack python
or python development on your system the installation fails. It seems to already be included in MAC OS X but
on other systems you should install/update it. Here is a 
<a href="http://wiki.python.org/moin/BeginnersGuide/Download" target="_blank" style="text-decoration:none;">guide.</a>
<br /> -->

Patchwork uses a python mod called pysam for "reading" the chromosomes. <br /> Currently 
<a href="https://code.google.com/p/pysam" target="_blank" style="text-decoration:none;">Pysam 0.7.4</a> 
and 
<a href="http://www.python.org/download/releases/2.6.7/" target="_blank" style="text-decoration:none;">Python 2.6.7 and 2.7</a>
are supported.
