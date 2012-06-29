karyotype_check <- function(chr,start,end,int,ai,Cn,mCn,t,name='',xlim=c(-1.02,1.02),ylim=0:1)
	{
	## TAPS scatter plot of a full sample, used for visual quality control. 
    png(paste(name,'_karyotype_check.png',sep=''),width=1300,height=1300)
 

    	#colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    	#colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")

    	#layout(matrix(1:25,nrow=5,byrow=T), widths=1,heights=1)
    
    	aix=!(is.na(mCn)|is.na(Cn))
    	chr=chr[aix]
    	start=start[aix]
    	end=end[aix]
    	int=int[aix]
    	ai=ai[aix]
    	Cn=Cn[aix]
    	mCn=mCn[aix]
    	pos <- (start+end)/2
    	length=end-start
    	labels=paste(Cn,mCn,sep='m')
    
    	size=rep(1,length(chr))
    	size[length>2000000]=2
    	size[length>5000000]=3
    	size[length>10000000]=4
    	
		## The variants present in this data
    	variants <- unique(labels)
    	variants_data <- NULL
    	
    	## Extracts the Log-R and Allelic Imbalance Ratio of variants
    	for (v in variants) 
    		{ 
        	t_cn <- paste('cn',strsplit(v,'m')[[1]][1],sep='')
        	t_int <- as.numeric(t$int[t_cn][[1]])
        	t_ai <- as.numeric(t$ai[paste('cn',v,sep='')][[1]])
        	variants_data <- rbind(variants_data,c(t_int,t_ai))
    		}
  
    	col <- rep('#B0B0B070',length(chr))
         
   	 	plot(int,ai,
        	pch=16,
        	cex=c(size),
        	cex.lab=3,
        	mar=c(0.1,0.1,0.1,0.1),
    		main = "",
			xlab = name,
			ylab = "",
			col = col,
			xlim = xlim,
			ylim = ylim,
			yaxt="n",
        	xaxt="n")
        	axis(1,cex.axis=2)
        	axis(2,cex.axis=2)
        	
		text(variants_data[,1],variants_data[,2],
        	labels=variants,
        	col='black',
        	cex=2)
    
    dev.off()
	}
