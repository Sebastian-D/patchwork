patchwork.createreference = function(...,output="REFOUT")
	{
	
	data=NULL
	
	files = list(...)
	
	data(ideogram,package="patchworkData")
	chroms=as.character(unique(ideogram$chr))
		
	j = 0
	
	for (i in files)
		{
		j = j + 1
		cat(paste("Reading Chromosomes of file ",i," \n",sep=""))
		data[[j]] = patchwork.readChroms(i,chroms)
		cat(paste("Initiating Gc Normalization of file ",i," \n",sep=""))
		data[[j]] = patchwork.GCNorm(data[[j]])
		data[[j]] = data[[j]][c("chr","pos",
											#"counts",
											"norm")]
		}
	
	if (j > 1)
		{
		cat("Merging the files into reference. \n")
		data_ = merge(data[[1]],data[[2]],all=TRUE,by=1:2)
		if (j > 2)
			{
			for (c in 3:j)
				{
				data_ = merge(data_,data[[c]],all=TRUE,by=1:2)
				}
			}
		#counts_ = seq(from=3,to=(length(data_) - 1),by=2)
		norm_ 	= seq(from=3,to=length(data_),by=1) 
	
		cat("Creating and applying mean functions to the reference. \n")
		normaldata = data.frame(chr=data_$chr,pos=data_$pos,
							#counts_mean=apply(data_[,counts_],
							#1,mean,na.rm=TRUE),
							reference=apply(data_[,norm_],
							1,mean,na.rm=TRUE))
		}
	else
		{
		data = data[[j]]
		cat("Creating and applying mean functions to the reference. \n")
		cat("Note: Creating from single bamfile. Cannot obtain mean ")
		cat("counts or norm however naming columns counts_mean and ")
		cat("norm_mean for compatability \n")
		normaldata = data.frame(chr=data$chr,pos=data$pos,#gc=data$gc,
							#gck=data$gck,gcm=data$gcm,
							#counts_mean=data$counts,
							reference=data$norm)
		}
	
	#conform chromosome naming to suit

	
	cat("Saving to ",output,".Rdata \n",sep="")
	save(normaldata,file=paste(output,".Rdata",sep=""))
	return(normaldata)
	}