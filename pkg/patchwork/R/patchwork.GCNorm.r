patchwork.GCNorm <- function(data)
	{
	#Load data included in package
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))
	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
	#load(paste(packagepath,"/data/normaldata.RData",sep=""))

	## GC content normalization
	info = data.frame(from=c(0,seq(10,90,2)),to=c(seq(10,90,2),100),mean=NA,median=NA,np=NA)

	#Remove NA's from gc and gck. I have only ever seen NA's in gck for HG18 genome.
	data = data[!is.na(data$gc),]
	data = data[!is.na(data$gck),]
	
	xnorm = data$counts
	
	#1st round
	for (i in 1:nrow(info))
		{
		ix = (data$gc >= info$from[i]) & (data$gc <= info$to[i])
		info$mean[i] = mean(data$counts[ix])
		info$median[i] = median(data$counts[ix])
		info$sd[i] = sd(data$counts[ix])
		info$np[i] = sum(ix)
		xnorm[ix] = 1 * ( data$counts[ix] / info$mean[i] )
		}	
		
	xnorm = round(xnorm,2)
		
	# 2nd round
	info2 = data.frame(from=c(0,seq(20,80,2)),to=c(seq(20,80,2),100),mean=NA,median=NA,np=NA)
	xnorm2 = xnorm
	
	for (i in 1:nrow(info2))
		{
		ix = (data$gck >= info2$from[i]) & (data$gck <= info2$to[i])
		info2$mean[i] = mean(xnorm[ix])
		info2$median[i] = median(xnorm[ix])
		info2$sd[i] = sd(xnorm[ix])
		info2$np[i] = sum(ix)
		xnorm2[ix] = 1 * ( xnorm[ix] / info2$mean[i] )
		}

	xnorm2 = round(xnorm2,2)
	data$norm=xnorm2
	return(data)
	}
