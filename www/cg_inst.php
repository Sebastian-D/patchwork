PatchworkCG works on Unix (Linux, Ubuntu, MacOSX, etc) based systems. You will probably get an error
when trying to install it on Windows, and you definitely can't run it.
<br /><br />

Begin by starting R. It is recommended that you use the latest version.

<br /><br />
Issue the following commands within R to install patchworkCG.

<pre>
    install.packages("patchworkCG", repos="http://R-Forge.R-project.org")
</pre>

If for some reason that does not work add the 'type="source"' to it, as so:

<pre>
    install.packages("patchworkCG", repos="http://R-Forge.R-project.org",type="source")
</pre>

If something still does not work with installation, here is a link to the location of the source files: <br />

<a href="http://r-forge.r-project.org/R/?group_id=1250" target="_blank" style="text-decoration:none">
	Patchwork</a>
