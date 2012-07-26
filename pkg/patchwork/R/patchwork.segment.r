patchwork.segment <- function(kbsegs,chroms,Alpha,SD)
	{
	#library(DNAcopy)
	#Load data included in package
	#packagepath = system.file(package="patchwork")
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))
	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
	#load(paste(packagepath,"/data/normaldata.RData",sep=""))
	data(ideogram,package="patchworkData")

	# Segmentation
	segs <- allsegs <- NULL;
	for (c in chroms)
		{
		psegments <- qsegments <- NULL
		ix= kbsegs$chr==c & kbsegs$pos < ideogram$start[ideogram$chr==c]
		if (sum(ix)>0)
			{ 
			psegments=segment(smooth.CNA(CNA(log2(kbsegs$ratio[ix]), c, kbsegs$pos[ix], data.type='logratio',sampleid=paste(c,'p')),smooth.region=40), undo.splits='sdundo', undo.SD=SD,min.width=5, alpha=Alpha)$output[,2:6]
			psegments$arm='p'
			}
		ix=kbsegs$chr==c & kbsegs$pos > ideogram$end[ideogram$chr==c]
		if (sum(ix)>0)
			{ 
			qsegments=segment(smooth.CNA(CNA(log2(kbsegs$ratio[ix]), c, kbsegs$pos[ix], data.type='logratio',sampleid=paste(c,'q')),smooth.region=40), undo.splits='sdundo', undo.SD=SD,min.width=5, alpha=Alpha)$output[,2:6]
			qsegments$arm='q'
			}
		segs=rbind(segs,psegments,qsegments)
		}	
	colnames(segs)=c('chr','start','end','np','mean','arm')
	return(segs)
	}
