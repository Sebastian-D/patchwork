patchwork.readChroms <- function(BamFile,chroms)
	{
	## Read chromosomes:
	# Ideogram data downloaded from UCSC table browser.
	
	#library(patchwork)
	#patchwork.plot("cancer.bam","pile.up",reference="illumina")
	
	#Load data included in package
	packagepath = system.file(package="patchwork")
	datapath = system.file(package="patchworkData")
	
	
	
	#if .python folder exists, delete it.
	test = getwd()
	try( setwd(".python"),silent=TRUE)
	if(getwd() != test)
		{
		setwd("..")
		system("rm -r .python")
		}
	
	
	#copy python folder
	system(paste("cp -r ",packagepath,"/python .python",sep=""))
	
	#If import pysam doesnt work, ie pysam isnt installed yet, we install pysam.
	if (system(paste("python ",getwd(),"/.python/Check.py",sep=""),intern=T) == "Pysam Not Available")
			{
			setwd(".python")
			cat("Pysam has not yet been installed on your system. \n Attempting to Build \n")
			system("sleep 5")
			# build pysam
			system("python setup.py build > pysam_build_output.txt")
			system("cp pysam_build_output.txt ..")
			cat("Build complete. Output of build can be viewed in file pysam_build_output.txt.
				\n Attempting to Install \n")
			# install pysam to site.USER_BASE
			system("python setup.py install --user > pysam_install_output.txt")
			system("cp pysam_install_output.txt ..")
			cat("Install complete. Output of install can be viewed in file pysam_install_output.txt. \n")
			setwd("..")
			}

	readlength=80 # this is currently hard coded and may be adjusted.
	seglength=200 # (do not change)
	data=NULL

	for (c in chroms) 
		{
		cat ('Reading',c,'\n')
		# Parse the bam file using python for each chromosome into a
		# new temporary file that holds position.
		system(paste("python .python/Pysamloader.py ",BamFile," ",c,sep =""))
		# Load the newly created position file
		unix = read.csv(".Tmp_chr_pos",sep='',header=F)$V1 + readlength/2

		#on chromosome 1, if the positions start below 10000 the hg18 has been
		#used for mapping. Otherwise its probably the hg19.
		if (c == chroms[1])
			{
			if (unix[1] <= 10000)
				{
				hg18 = TRUE
				}
			else
				{
				hg18 = FALSE
				}
			}

		if (hg18==TRUE)
			{
			gc = readRDS(paste(datapath,"/extdata/hg18.",c,".mark.rds",sep=""))
			}
		else
			{
			gc = readRDS(paste(datapath,"/extdata/",c,".mark.rds",sep=""))
			}

		gc = gc[!is.na(gc$pos) & !is.na(gc$gc250),]
		gc = gc[order(gc$pos),]
		xpos=gc[,1]
		xgc=gc[,2]
		xgck=gc[,3]
		xgcm=gc[,4]
		xn=length(xpos)
		xbreaks=c(xpos[1]-seglength,xpos)+seglength/2
		xseglength=xbreaks[2:length(xbreaks)]-xbreaks[1:(length(xbreaks)-1)]
		xhist=hist(unix[unix>xbreaks[1]&unix<xbreaks[length(xbreaks)]],breaks=xbreaks,plot=F)
	
		ix=xseglength==seglength & xhist$counts>0

		xpos=xhist$mids[ix]
		xcounts=xhist$counts[ix]
		xgc=xgc[ix]
		xgck=xgck[ix]
		xgcm=xgcm[ix]
		xn=sum(ix)
		data=rbind(data,data.frame(chr=c,pos=xpos,gc=xgc,gck=xgck,gcm=xgcm,counts=xcounts))
		}
	system("rm -r .python")
	system("rm .Tmp_chr_pos")
	return(data)
	}
