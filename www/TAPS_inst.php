TAPS works on Unix (Linux, MacOSX, etc) based systems and has been tested with Rstudio on Windows.
<br /><br />

Begin by starting R. It is recommended that you use the latest version.
<br /><br />
Enable all repositories in R so the dependencies of TAPS can be found.<br /><br />

<pre>
    setRepositories()
</pre>
Install the packages that TAPS depends on to function correctly. <br /><br />

<pre>
    source("https://bioconductor.org/biocLite.R")
    biocLite("DNAcopy")
    biocLite("affxparser")
    install.packages("fields")
    install.packages("foreach")
    install.packages("jpeg")
    install.packages("xlsx")

    #Only on unix-based systems
    install.packages("doMC")
</pre>

Install TAPS.

<br /><br />

<pre>
    install.packages("TAPS", repos="http://R-Forge.R-project.org")
</pre>

If that does not work, try installing it from source by adding 'type="source"' to the command.
<pre>
    install.packages("TAPS", repos="http://R-Forge.R-project.org",type="source")
</pre>

If something still does not work with installation, here is a link to the location of the source file: <br />

<a href="http://r-forge.r-project.org/R/?group_id=1250" target="_blank" style="text-decoration:none;">
	TAPS</a><br />
