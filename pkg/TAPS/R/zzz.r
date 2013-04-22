.onLoad <- function(...){
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to TAPS")
	packageStartupMessage("\n")
	packageStartupMessage("Version: ",utils::packageDescription('TAPS')$Version)
	packageStartupMessage("\n")
	packageStartupMessage("If this is your first time running TAPS you should visit the")
	packageStartupMessage("homepage at (http://patchwork.r-forge.r-project.org/) or see ?TAPS.")
	packageStartupMessage("\n")
	packageStartupMessage("Remember to keep TAPS up to date: update.packages(repos=\"http://R-Forge.R-project.org\")")
}