CG_KaCh <- function(chr,start,end,int,ai,mchr,mpos,mval,
					sval,schr,spos,name='',xlim=c(0,2.4),ylim=c(0.3,1))  
	{   
    
    #packagepath = system.file(package="patchwork")
 	
    #load(paste(packagepath,"/data/ideogram.RData",sep=""))
    
    data(ideogram,package="patchworkCG")
    
    #create a color ramp palette to segue the color scheme for plot.
    
    #colors_p <- colorRampPalette(c("#2754e8","#ad09e8"),space="rgb")
    #colors_q <- colorRampPalette(c("#c309e8","#e81e2c"),space="rgb")

    colors_p <- colorRampPalette(c("#2754e8","#ad09e8"),space="rgb")
    colors_q <- colorRampPalette(c("#c309e8","#e81e2c"),space="rgb")

    ai[is.na(ai)]=0
    aix=ai!=0
    chr=chr[aix]
    start=start[aix]
    end=end[aix]
    int=int[aix]
    ai=ai[aix]
    pos <- (start+end)/2
    length=end-start
    
    size=rep(1,length(chr))
    size[length>1000000]=2
    size[length>3000000]=3
    size[length>5000000]=4
    size[length>10000000]=5
    size[length>20000000]=6

    for (c in 1:24) 
    	{
        this <- ideogram[ideogram$c==c,]
        
        png(paste(name,'_KaCh_',this$chr,'.png',sep=''),width=1080,height=1300)
        	layout(matrix(1:3,nrow=3,byrow=T), widths=10,heights=c(10,5,5))
        	
        	mar.default <- c(5,4,4,2) + 0.1
			par(mar = mar.default + c(0, 4, 0, 0)) 
        
       	 	ix <- chr==this$chr
            ix_p = ix & (pos < this$start)
            ix_q = ix & (pos > this$end)
            col <- rep('#B0B0B070',length(chr))

            #loc = list of colors
            loc_p=colors_p(ceiling(this$start/100))
            loc_q=colors_q(ceiling((this$length - this$end)/100))

            col[ix_p] = paste(loc_p[ceiling(pos[ix_p]/100)],'70',sep='')
            col[ix_q] = paste(loc_q[ceiling((pos[ix_q]-this$end)/100)],'70',sep='')
                   
        	#col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')  
        	#col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='') 
        
        	plot(c(int[!ix],int[ix]),c(ai[!ix],ai[ix]),
            	pch=16,
            	cex=c(size[!ix],size[ix]),
            	cex.lab=3,
            	#mar=c(0.1,0.1,0.1,0.1),
		    	main = "",
		    	xlab = "",
		    	ylab = '',
		    	col = c(col[!ix],col[ix]),
		    	xlim = xlim,
		    	ylim = ylim,
		    	yaxt="n",
        		xaxt="n")
        	axis(1,at=seq(0,2.5,0.1),cex.axis=2)
        	axis(2,at=seq(0,1,0.1),cex.axis=2)
        	mtext(text="Allelic imbalance",side=2,line=4,las=3,cex=2)
        	mtext(text="Coverage, all segments",side=1,line=4,cex=2,las=1)
		
        	par(new=F)
        	mix <- mchr==as.character(this$chr)
        	#col <- rep('#D0D0D0D0',sum(mix))       
        	#col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')  
        	#col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='') 
        	col=rep('#000000',sum(ix))
            #Remove Alpha.
            #col[ix_p] = loc_p[ceiling(pos[ix_p]/100)]
            #col[ix_q] = loc_q[ceiling((pos[ix_q]-this$end)/100)]
            col[pos[ix] < this$mid] <- loc_p[ceiling(pos[ix][pos[ix] < this$mid]/100)]
            col[pos[ix] > this$mid] <- loc_q[ceiling((pos[ix][pos[ix] > this$mid] - this$mid)/100)]

            #Origional
        	#col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
        	#col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))
        
        	plot(mpos[mix],mval[mix],
            	pch=16,
            	cex=2,
            	cex.lab=3,
            	#mar=c(0.1,0.1,0.1,0.1),
    	    	main = "",
		    	xlab = "",
		    	ylab = "",
		    	col = '#00000010',
		    	xlim = c(0,this$length),
		    	ylim = xlim,
		    	yaxt="n",
        		xaxt="n")
        	axis(1,cex.axis=2)
        	axis(2,cex.axis=2)
        	mtext(text=paste("Coverage,",this$chr),side=1,line=4,cex=2,las=1)
		    
        	segments(x0=start[ix],x1=end[ix],
        	    y0=int[ix],y1=int[ix],                
            	col=col,
            	lwd=4,)
            
       	 	par(new=F)       
        	six <- schr==as.character(this$chr)
        	plot(spos[six],sval[six],
            	pch=16,
            	cex=2,
            	cex.lab=3,
            	#mar=c(0.1,0.1,0.1,0.1),
            	main = "",
		    	xlab = "",
		    	ylab = "",
		    	col = '#00000010',
		    	xlim = c(0,this$length),
		    	ylim = c(0,1),
		    	yaxt="n",
        		xaxt="n")
        	axis(1,cex.axis=2)
        	axis(2,cex.axis=2)
        	mtext(text=paste("Allelic imbalance,",this$chr),side=1,line=3.5,cex=2,las=1)
        	
        	segments(x0=start[ix],x1=end[ix],
        	    y0=ai[ix],y1=ai[ix],                
            	col=col,
            	lwd=4,)
            	
    	dev.off()
    	}
	}
