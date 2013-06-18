TAPS works on Unix (Linux, Ubuntu, MacOSX, etc) based systems and has been tested with Rstudio on Windows.
<br /><br />

Begin by starting R. It is recommended that you use the latest version.
<br /><br />
Start by enabling all repositories in R so the dependencies of TAPS can be found.<br /><br />

<pre>
    setRepositories()
</pre>
Install the packages that TAPS depends on to function correctly. <br /><br />

<pre>
    install.packages("DNAcopy")
    install.packages("affxparser")
    install.packages("fields")
</pre>

To install TAPS enter this command in R.

<br /><br />

<pre>
    install.packages("TAPS", repos="http://R-Forge.R-project.org")
</pre>

If for some reason that does not work try installing it from source by adding 'type="source"' to the command.
<pre>
    install.packages("TAPS", repos="http://R-Forge.R-project.org",type="source")
</pre>

If something still does not work with installation, here are links to the location of the source files: <br />

<a href="http://r-forge.r-project.org/R/?group_id=1250" target="_blank" style="text-decoration:none;">
	TAPS</a><br />