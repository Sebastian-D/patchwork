.onLoad <- function(...){
	#supressMessages(library(DNAcopy))
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to patchwork.")
	packageStartupMessage("\n")
	packageStartupMessage("Version: ",utils::packageDescription('patchwork')$Version)
	packageStartupMessage("\n")
	packageStartupMessage("If this is your first time running patchwork you should visit the")
	packageStartupMessage("homepage at (http://patchwork.r-forge.r-project.org/) or see ?patchwork.readme")
	packageStartupMessage("or the README file found at the install location of patchwork/inst/README")
	}
