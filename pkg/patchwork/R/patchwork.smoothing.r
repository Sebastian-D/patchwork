patchwork.smoothing <- function(data,normaldata=NULL,reference=NULL,chroms)
	{
	#Load data included in package
	#packagepath = system.file(package="patchwork")
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))
	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
	#load(paste(packagepath,"/data/normaldata.RData",sep=""))
	data(ideogram,package="patchworkData")	

	if (is.null(normaldata)) load(reference)	

	data <- merge(data,normaldata,all=T)
	data <- data[order(data$pos),]; data <- data[order(data$chr),]
	colnames(data)[8] <- 'reference'

	# apply smoothing @ 10kb (~50 minibins)
	seglength <- 10000
	window <- seglength/200 
	kbsegs=NULL
	for (c in chroms) 
		{
		cat("Smoothing Chromosome:",c,'\n')

		# p-arm
		i=which(as.character(ideogram$chr)==c) #!!
		
		#attach(data[as.character(data$chr)==c & data$pos<ideogram$start[i],])
		
		ix = as.character(data$chr)==c & data$pos<ideogram$start[i]
		norm = data$norm[ix]
		pos = data$pos[ix]
		reference = data$reference[ix]
		
		n=length(pos)
		if (n>1000) 
			{ 
			# p-arm may be missing.
			kbposix <- seq(window,n-window,by=window)
			kbpos <- pos[kbposix]
			nsegs <- length(kbpos)
			coverage <- rep(NA,nsegs)
			refcoverage <- rep(NA,nsegs)

			for (i in 1:nsegs) 
				{
				j=kbposix[i]
				ix <- (j-(window/2)):(j+(window/2))
				coverage[i] <- median(norm[ix], na.rm=T)
				refcoverage[i] <- median(reference[ix], na.rm=T)
				}
			kbsegs <- rbind(kbsegs,data.frame('chr'=c,kbpos,coverage,refcoverage))
			}
			
		#detach()
	
		# q-arm
		i=which(as.character(ideogram$chr)==c) #!!
		#attach(data[as.character(data$chr)==c & data$pos>ideogram$end[i],])
		
		ix = as.character(data$chr)==c & data$pos>ideogram$end[i]
		norm = data$norm[ix]
		pos = data$pos[ix]
		reference = data$reference[ix]
		
		n=length(pos)
		kbposix <- seq(window,n-window,by=window)
		kbpos <- pos[kbposix]
		nsegs <- length(kbpos)
		coverage <- rep(NA,nsegs)
		refcoverage <- rep(NA,nsegs)

		for (i in 1:nsegs) 
			{
			j=kbposix[i]
			ix <- (j-(window/2)):(j+(window/2))
			coverage[i] <- median(norm[ix], na.rm=T)
			refcoverage[i] <- median(reference[ix], na.rm=T)
			}
			
		kbsegs <- rbind(kbsegs,data.frame('chr'=c,kbpos,coverage,refcoverage))
		#detach()
		}
	colnames(kbsegs) <- c('chr','pos','coverage','refcoverage')
	kbsegs$ratio <- kbsegs$coverage/kbsegs$refcoverage
	return(kbsegs)
	}
