.onLoad <- function(...){
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to patchwork for CompleteGenomics data!")
	packageStartupMessage("\n")
	packageStartupMessage("Version: ",utils::packageDescription('patchworkCG')$Version)
	packageStartupMessage("\n")
	packageStartupMessage("If this is your first time, please see the documentation at the")
	#packageStartupMessage("the README file at the default installation location, the")
	packageStartupMessage("homepage (http://patchwork.r-forge.r-project.org/) or in ?patchworkCG.readme")
	packageStartupMessage("or the README file found at the install location of patchworkCG/example_plots/README")
	}