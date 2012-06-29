CG_karyotype <- function(chr,start,end,int,ai,
						name='',xlim=c(0,2.4),ylim=c(0.3,1)) 
	{   
    png(paste(name,'_karyotype.png',sep=''),width=1300,height=1300)
 
    	#packagepath = system.file(package="patchwork")
 	
    	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
    	
    	data(ideogram,package="patchworkCG")
    	
		#colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    	#colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")

        colors_p <- colorRampPalette(c("#2754e8","#ad09e8"),space="rgb")
        colors_q <- colorRampPalette(c("#c309e8","#e81e2c"),space="rgb")

    	layout(matrix(1:25,nrow=5,byrow=T), widths=1,heights=1)
    	
    	aix=ai!=0
    	chr=chr[aix]
    	start=start[aix]
    	end=end[aix]
    	int=int[aix]
    	ai=ai[aix]
    	pos <- (start+end)/2
    	length=end-start
    
    size=rep(0.5,length(chr))
    size[length>1000000]=1
    size[length>3000000]=2
    size[length>5000000]=3
    size[length>10000000]=4
    size[length>20000000]=5
 
    
    	for (c in 1:24) 
    		{

        	this <- ideogram[ideogram$c==c,]
        	ix <- chr==this$chr
        	ix[is.na(ix)]=FALSE

            ix_p = ix & (pos < this$start)
            ix_q = ix & (pos > this$end)
            col <- rep('#B0B0B070',length(chr))

            #loc = list of colors
            loc_p=colors_p(ceiling(this$start/100))
            loc_q=colors_q(ceiling((this$length - this$end)/100))

            col[ix_p] = paste(loc_p[ceiling(pos[ix_p]/100)],'70',sep='')
            col[ix_q] = paste(loc_q[ceiling((pos[ix_q]-this$end)/100)],'70',sep='')
 
        	plot(c(int[!ix],int[ix]),c(ai[!ix],ai[ix]),
        		pch=16,
        		cex=c(size[!ix],size[ix]),
        		cex.lab=3,
        		mar=c(0.1,0.1,0.1,0.1),
				main = "",
				xlab = this$chr,
				ylab = "",
				col = c(col[!ix],col[ix]),
				xlim = xlim,
				ylim = ylim,
		    	yaxt="n",
        		xaxt="n")
        	axis(1,cex.axis=1.5)
        	axis(2,cex.axis=1.5)
			par(new=F)  
    		}
    dev.off()
	}


    
    



