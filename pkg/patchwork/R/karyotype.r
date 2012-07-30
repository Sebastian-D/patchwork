karyotype <- function(chr,start,end,int,ai,
						name='',xlim=c(0,2.4),ylim=c(0,1)) 
	{   
    png(paste(name,'.karyotype.png',sep=''),width=1300,height=1300)
 
    	#packagepath = system.file(package="patchwork")
 	
    	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
    	
    	data(ideogram,package="patchworkData")
    	
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
    	size[length>2000000]=1
    	size[length>5000000]=2
    	size[length>10000000]=3
 
    
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

            #In some cases, it seems HG18, we get an out of bounds for the segments colors.
            #This results in a NA which calls the color NA70 which obviously does not exists.
            #To solve this i catch these extremes and simply give them the most extreme color.
            #It should only possibly happen at the ends as loc_Q/loc_p start at 1.

            #End of loc_p failed
            p = length(col[ix_p])
            if(p >= 1)
                {
                while(col[ix_p][p]=="NA70")
                    {
                    col[ix_p][p] = paste(loc_p[length(loc_p)],'70',sep='')
                    p = p - 1
                    }
                }

            #End of loc_q failed
            q = length(col[ix_q])
            if(q >= 1)
                {
                while(col[ix_q][q]=="NA70")
                    {
                    col[ix_q][q] = paste(loc_q[length(loc_q)],'70',sep='')
                    q = q - 1
                    }
                }

 
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


    
    



