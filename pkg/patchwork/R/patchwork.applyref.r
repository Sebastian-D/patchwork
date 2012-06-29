patchwork.applyref <- function(data,reference=NULL)
	{

	if (ncol(data)>7) data <- data[,1:7]

	data <- merge(data,reference,all=F,by=1:2) ### Should change to all=T, but that must be fixed in smoothing.
	colnames(data)[8] = 'reference'
	data=data[order(data$pos),]
	data=data[order(data$chr),]
	return(data)
	}