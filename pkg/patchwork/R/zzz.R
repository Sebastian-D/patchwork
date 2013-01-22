.onLoad <- function(...){
	#supressMessages(library(DNAcopy))
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to patchwork.")
	packageStartupMessage("\n")
	# This was an attempt to check if patchwork needed updating. It failed. Running python in onLoad is not trivial.
	if(system(paste("python ",system.file(package="patchwork"),"/python/vercheck.py",sep=""),intern=T) != utils::packageDescription('patchwork')$Version)
		{
		packageStartupMessage("Patchwork is not up to date, please update at:")
		packageStartupMessage("http://patchwork.r-forge.r-project.org/")
	  packageStartupMessage("or using: update.packages(repos=\"http://R-Forge.R-project.org\")")
		packageStartupMessage("\n")
		}
	else
		{
		packageStartupMessage("Version: ",utils::packageDescription('patchwork')$Version)
		packageStartupMessage("\n")
		packageStartupMessage("If this is your first time running patchwork you should visit the")
		packageStartupMessage("homepage at (http://patchwork.r-forge.r-project.org/) or see ?patchwork.readme")
		packageStartupMessage("or the README file found at the install location of patchwork/inst/README")
		packageStartupMessage("\n")
		packageStartupMessage("Remember to keep patchwork up to date: update.packages(repos=\"http://R-Forge.R-project.org\")")
		}
	}
