# karyotype_chromsCN <- function(chr,start,end,int,ai,Cn,mCn,mchr,mpos,mval,
# 								schr,spos,sval,name='',xlim=c(-1.02,1.82),ylim=0:1,
# 								maxCn=8)  
# 	{   
    
#     #packagepath = system.file(package="patchwork")
 	
#     #load(paste(packagepath,"/data/ideogram.RData",sep=""))
    
#     data(ideogram,package="patchworkData")
    
#     colors_p <- colorRampPalette(c("#2754e8","#ad09e8"),space="rgb")
#     colors_q <- colorRampPalette(c("#c309e8","#e81e2c"),space="rgb")

#     ai[is.na(ai)]=0
#     aix=ai!=0
#     chr=chr[aix]
#     start=start[aix]
#     end=end[aix]
#     int=int[aix]
#     ai=ai[aix]
#     Cn=Cn[aix]
#     mCn=mCn[aix]
#     pos <- (start+end)/2
#     length=end-start
    
#     size=rep(1,length(chr))
#     size[length>2000000]=2
#     size[length>5000000]=3
#     size[length>10000000]=4
       
#     for (c in 1:24) 
#     	{
#         this <- ideogram[ideogram$c==c,]
        
#         png(paste(name,'_karyotypeCN.',this$chr,'.png',sep=''),width=1080,height=1300)
#         	layout(matrix(1:4,nrow=4,byrow=T), widths=10,heights=c(10,5,5,5))
        	
#         	mar.default <- c(5,4,4,2) + 0.1
# 			par(mar = mar.default + c(0, 4, 0, 0)) 
        	
#             ix <- chr==this$chr
#             ix_p = ix & (pos < this$start)
#             ix_q = ix & (pos > this$end)
#             col <- rep('#B0B0B070',length(chr))

#             #loc = list of colors
#             loc_p=colors_p(ceiling(this$start/100))
#             loc_q=colors_q(ceiling((this$length - this$end)/100))

#             col[ix_p] = paste(loc_p[ceiling(pos[ix_p]/100)],'70',sep='')
#             col[ix_q] = paste(loc_q[ceiling((pos[ix_q]-this$end)/100)],'70',sep='')

#             #In some cases, it seems HG18, we get an out of bounds for the segments colors.
#             #This results in a NA which calls the color NA70 which obviously does not exists.
#             #To solve this i catch these extremes and simply give them the most extreme color.
#             #It should only possibly happen at the ends as loc_Q/loc_p start at 1.

#             #End of loc_p failed
#             p = length(col[ix_p])
#             if(p >= 1)
#                 {
#                 while(col[ix_p][p]=="NA70")
#                     {
#                     col[ix_p][p] = paste(loc_p[length(loc_p)],'70',sep='')
#                     p = p - 1
#                     }
#                 }

#             #End of loc_q failed
#             q = length(col[ix_q])
#             if(q >= 1)
#                 {
#                 while(col[ix_q][q]=="NA70")
#                     {
#                     col[ix_q][q] = paste(loc_q[length(loc_q)],'70',sep='')
#                     q = q - 1
#                     }
#                 }
        
#         	plot(c(int[!ix],int[ix]),c(ai[!ix],ai[ix]),
#             	pch=16,
#             	cex=c(size[!ix],size[ix]),
#             	cex.lab=3,
#             	#mar=c(0.1,0.1,0.1,0.1),
#     	    	main = "",
# 		    	xlab = "",
# 		    	ylab = '',
# 		    	col = c(col[!ix],col[ix]),
# 		    	xlim = xlim,
# 		    	ylim = ylim,
# 		    	yaxt="n",
#         		xaxt="n")
#         	axis(1,at=seq(0,2.5,0.1),cex.axis=2)
#         	axis(2,at=seq(0,1,0.1),cex.axis=2)
#         	mtext(text="Allelic imbalance",side=2,line=4,las=3,cex=2)
#         	mtext(text="Coverage, all segments",side=1,line=4,cex=2,las=1)
#             if(c==1)
#                 {
#                 title(main=paste("Sample: ",name,sep=""),cex.main=2)
#                 }
		    	
# 			col=rep('#000000',sum(ix))
#         	col[pos[ix] < this$mid] <- loc_p[ceiling(pos[ix][pos[ix] < this$mid]/100)]
#             col[pos[ix] > this$mid] <- loc_q[ceiling((pos[ix][pos[ix] > this$mid] - this$mid)/100)]
#         	par(new=F)
        	
#         	plot(1,1,type='n',
#             	xlab="",
#             	ylab='',
#             	cex.lab=3,
#             	xlim = c(0,this$length),
#             	ylim = c(-0.1,8.1),
# 		    	yaxt="n",
#         		xaxt="n")
#         	axis(1,cex.axis=2)
#         	axis(2,cex.axis=2)
#         	mtext(text=paste('Total and minor copy number', this$chr),side=1,line=4,cex=2,las=1)
            	
#         	segments(x0=start[ix],x1=end[ix],
#          	   y0=Cn[ix],y1=Cn[ix],                
#             	col=col,
#             	lwd=5,)   #Om fel, prova ta bort detta komma
            	
#         	segments(x0=start[ix],x1=end[ix],
#             	y0=mCn[ix],y1=mCn[ix],                
#             	col='#BBBBBB',
#             	lwd=5,) #Om fel, prova ta bort detta komma
 
#         	par(new=F)
#         	mix <- as.character(mchr)==this$chr        
#         	plot(mpos[mix],mval[mix],
#             	pch=16,
#             	cex=2,
#             	cex.lab=3,
#             	#mar=c(0.1,0.1,0.1,0.1),
#     	    	main = "",
# 		    	xlab = "",
# 		    	ylab = "",
# 		    	col = '#00000010',
# 		    	xlim = c(0,this$length),
# 		    	ylim = xlim,
# 		    	yaxt="n",
#         		xaxt="n")
#         	axis(1,cex.axis=2)
#         	axis(2,cex.axis=2)
#         	mtext(text=paste("Coverage,",this$chr),side=1,line=4,cex=2,las=1)
		    	
#         	segments(x0=start[ix],x1=end[ix],
#             	y0=int[ix],y1=int[ix],                
#             	col=col,
#             	lwd=4,)
            	
#         	par(new=F)       
#         	six <- as.character(schr)==this$chr
        	
#         	plot(spos[six],sval[six],
#         	    pch=16,
#             	cex=2,
#             	cex.lab=3,
#             	#mar=c(0.1,0.1,0.1,0.1),
#             	main = "",
# 		    	xlab = "",
# 		    	ylab = "",
# 		    	col = '#00000010',
# 		    	xlim = c(0,this$length),
# 		    	ylim = c(0,1),
# 		    	yaxt="n",
#         		xaxt="n")
#         	axis(1,cex.axis=2)
#         	axis(2,cex.axis=2)
#         	mtext(text=paste("Allelic imbalance,",this$chr),side=1,line=3.5,cex=2,las=1)

#             segments(x0=start[ix],x1=end[ix],
#                 y0=ai[ix],y1=ai[ix],                
#                 col=col,
#                 lwd=4,)
#         dev.off()
#     	}
# 	}

#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#


karyotype_chromsCN <- function(chr,start,end,int,ai,Cn,mCn,
    mchr,mpos,mval,
    schr,spos,sval,t,
    name='',parameters,xlim=c(0,2.5),ylim=0:1,maxCn=8)  
{   
    
    #Get ideogram
    #data(ideogram,package="patchworkData")
    
    #Set color gradient for p and q arm of chromosome
    colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")
    
    #filter the data (remove NA) and set some standard parameters
    ai[is.na(ai)]=0
    aix=ai!=0
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
    
    #Set sizes for circles in whole genome plot
    size=rep(0.5,length(chr))
    size[length>2000000]=1
    size[length>5000000]=1.5
    size[length>10000000]=2
    
    # The variants present in this data
    variants <- unique(labels)
    variants_data <- NULL
    # Extracts the Log-R and Allelic Imbalance Ratio of variants    
    for (v in variants) 
    {
        t_cn <- paste('cn',strsplit(v,'m')[[1]][1],sep='')
        t_int <- as.numeric(t$int[t_cn][[1]])
        t_ai <- as.numeric(t$ai[paste('cn',v,sep='')][[1]])
        variants_data <- rbind(variants_data,c(t_int,t_ai))
    }
    
    #Loop over each chromosome. (Not chromosomeY as of writing)     
    for (c in 1:23) 
    {
        #Select the chromosome
        this <- ideogram[ideogram$c==c,]
        
        #Initialize jpeg
        jpeg(paste(name,'_karyotypeCN.',this$chr,'.jpg',sep=''),width=11.7,height=8.3,units="in",res=300)
        
        #split plot into desired formation
        split.screen(as.matrix(data.frame(left=c(0.05,rep(0.47,4)),
                                          right=c(0.45,rep(1,4)),
                                          bottom=c(0.25,0.63,0.42,0.37,0.17),
                                          top=c(0.75,0.83,0.62,0.42,0.37)))) #screen 1-5
        
        
        #Screen configuration overview
        #-------------------------------------------------------------
        #-------------------------------------------------------------
        #                             |                              |
        #                             |              2               |
        #                             |                              |
        #                             |------------------------------| 
        #                             |                              |
        #                             |              3               |
        #              1              |                              |
        #                             |------------------------------| 
        #                             |              4               |
        #                             |------------------------------|
        #                             |                              |
        #                             |              5               |
        #                             |                              |
        #-------------------------------------------------------------
        
        
        #------------------------------------------------------------
        #Left side of the plot (Whole genome)
        #------------------------------------------------------------
        #index of the selected chromsome to the chr object
        ix <- chr==as.character(this$chr)
        
        notix1 = chr=="chrX" & int > 2.5
        notix = !ix + notix1
        
        #Create an index of colors relating to the positions and lengths on this chromosome
        col <- rep('#46464620',length(chr))       
        col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')  
        col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='')
        
        #Go to screen 1
        screen(1)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        if(c!=23)
        {
            #Whole genome overview plot
            #Note that we are plotting chrY in the background, despite the fact that it is not an "active" chromosome.
            plot(c(int[notix],int[ix]),c(ai[notix],ai[ix]),
                 pch=16,
                 cex=c(size[notix],size[ix]),
                 main = "",
                 xlab = "",
                 ylab = "",
                 axes=F, #Remove axis
                 col = c(col[notix],col[ix]),
                 xlim = c(0,2.5),
                 ylim = ylim)
        }
        else
        {
            plot(int[notix],ai[notix],
                 pch=16,
                 cex=size[notix],
                 main = "",
                 xlab = "",
                 ylab = "",
                 axes=F, #Remove axis
                 col = col[notix],
                 xlim = c(0,2.5),
                 ylim = ylim)
        }
        
        xix=chr=="chrX"
        if(c!=23)
        {
            points(int[xix],ai[xix],pch=4,col='#46464680')
        }
        else
        {
            points(int[xix],ai[xix],pch=4,col=col[ix]) 
        }
        
        #Insert Y and X axis
        axis(side=2,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),las=1)
        axis(side=1,cex.axis=0.6,tck=0.926,col.ticks='#808080',at=seq(from=0,to=2.5,by=0.1),pos=0)
        
        #Add cn and mcn labels
        text(variants_data[,1],variants_data[,2],
             labels=variants,
             col='black',
             cex=0.75)
        
        #Titles and axis labels
        mtext(text="Normalized coverage",side=1,line=0.6,cex=1)
        mtext(text="Allelic imbalance",side=2,line=1,cex=1)
        
        #Check if name of sample is too long to be graphically pleasing and cut it if that is the case
        if(nchar(name)>12)
        {
            mtext(paste("Detailed view with copynumbers of sample: ",substring(name,1,12),"\n Chromsome ",c,sep=""),side=3)
        }
        else
        {
            mtext(paste("Detailed view with copynumbers of sample: ",name,"\n Chromsome ",c,sep=""),side=3)
        }
        
        #------------------------------------------------------------
        #Right side, Top (Total and minor copy numbers)
        #------------------------------------------------------------
        
        #Go to screen 2
        screen(2)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        #Set colors for total and minor copynumber
        col=rep('#000000',sum(ix))
        col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
        col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))
        
        #Empty plot for total and minor copynumbers
        plot(1,1,type='n',
             xlab='',
             ylab='',
             axes=F,
             xlim = c(0,this$length),
             ylim = c(-0.1,8.1))
        
        #Add total copnumber segments
        segments(x0=start[ix],x1=end[ix],
                 y0=Cn[ix],y1=Cn[ix],                
                 col=col,
                 lwd=5)
        
        #Add minor copynumber segments
        segments(x0=start[ix],x1=end[ix],
                 y0=mCn[ix],y1=mCn[ix],                
                 col='#BBBBBB',
                 lwd=5)
        
        #Add Y axis
        axis(side=2,tck=0.926,col.ticks='#808080',at=seq(0,8,by=1),cex.axis=0.6,pos=0,las=1)
        axis(side=4,labels=F,tck=0,pos=this$length)
        
        #Add Y label and date/time title
        mtext("Total and minor copynumber",side=2,line=0.3)
        mtext(format(Sys.time(),"%Y-%m-%d %H:%M"),side=3,cex=0.6,adj=0.95)
        
        #------------------------------------------------------------
        #Right side, Middle Top - Signal
        #------------------------------------------------------------
        #Go to screen 3
        screen(3)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        #Select the correct chromosome and remove stuff lower than -1 
        mix <- mchr==as.character(this$chr) #& mval>(-1)
        
        #Predefine ymin ymax and sequence between them
        # ymin=floor(min(int[ix]))-0.5
        # if(ymin > -1)
        # {
        #     ymin = -1
        # }
        # ymax=ceiling(max(int[ix]))+0.5
        # if(ymax < 1)
        # {
        #     ymax = 1
        # }
        ymin = 0
        ymax = 2.5
        seqminmax=seq(ymin,ymax,by=0.5)
        
        #Create an index of colors relating to the positions and lengths on this chromosome
        col=rep('#000000',sum(ix))
        col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
        col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))
        
        #Plot log-ratio over position
        plot(mpos[mix],mval[mix],
             pch=20,
             cex=0.5,
             main = "",
             xlab = "",
             ylab = "",
             axes=F,
             col = '#00000005',
             xlim = c(0,this$length),
             ylim = c(ymin,ymax))
        
        #Add colored segments based on the log-ratio data
        segments(x0=start[ix],x1=end[ix],
                 y0=int[ix],y1=int[ix],                
                 col=col,
                 lwd=4)
        
        #Add X and Y axis. The (side=4) axis is just a black line showing where the data ends
        axis(side=2,tck=0.926,col.ticks='#808080',at=seqminmax, #seq(from=-1,to=1.5,by=0.5),
             #labels=c("-1","-.5","0",".5","1","1.5"),
             cex.axis=0.6,pos=0,las=1)
        axis(side=4,labels=F,tck=0,pos=this$length)
        axis(side=1,tck=0.926,col.ticks='#808080',at=seq(5e6,this$length,by=5e6),
             labels=FALSE,cex.axis=0.6,pos=ymin)
        
        #Add Y axis label
        mtext("Norm. coverage",side=2,line=0.3)
        
        #------------------------------------------------------------
        #Right side, Middle Bottom - Cytobands
        #------------------------------------------------------------
        
        #Go to screen 4
        screen(4)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
        plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,0.3),xlim=c(0,max(mpos[mix])))
        
        #Add Y axis label
        mtext(text="Cytoband",side=2,las=1,line=-1,cex=0.8)
        
        #index of the chromData for this chromosome
        dix = chromData$chr == as.character(this$chr)
        
        #Add cytoband information as differently colored rectangles
        rect(xleft=chromData$chromStart[dix],xright=chromData$chromEnd[dix],
             #ybottom=0.15-chromData$thickness[dix]*0.03,ytop=0.15+chromData$thickness[dix]*0.03,
             ybottom=0,ytop=0.3,
             col=paste(chromData$col[dix],'99',sep=''),
             lty=0)
        
        #Att text annotation to cytobands
        for(j in 1:length(chromData$name[dix]))
        {
            #If the cytoband region is larger than 5 Mb, it can fit some text
            if((chromData$chromEnd[dix][j] - chromData$chromStart[dix][j]) > 4*10^6)
            {
                # Text is added at x position (all of previous chromosome) + (middle of region)
                # The "srt" flips the text 90 degrees
                text(x=(chromData$chromStart[dix][j]+((chromData$chromEnd[dix][j]-chromData$chromStart[dix][j])/2)),
                     y=0.15,labels=chromData$name[dix][j], srt=90,cex = 0.5,xpd=T)
            }
        }
        
        
        #------------------------------------------------------------
        #Right side, Bottom - AI
        #------------------------------------------------------------ 
        
        #Go to screen 5
        screen(5)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        #Index of correct chromsome and allele frequency not in either 0 or 1
        six <- schr==as.character(this$chr) & !(sval %in% c(0,1))
        
        #plot allele frequency over position
        plot(spos[six],sval[six],
             pch=20,
             cex=0.5,
             main = "",
             xlab = "",
             ylab = "",
             #yaxt="n",
             axes=F,
             col = '#00000010',
             xlim = c(0,this$length),
             ylim = c(0,1))
        
        #Add X and Y axis as well as a line to the rightmost of the data (side=4)
        axis(side=2,tck=0.926,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=0,las=1)
        axis(side=4,labels=F,tck=0,pos=this$length)
        axis(side=1,tck=0.925,col.ticks='#808080',at=seq(5e6,this$length,by=5e6),
             labels=paste(seq(5,round(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
        
        #Add X and Y label
        mtext("Allelic imbalance",side=2,line=0.3)
        mtext("Position (Mb)",side=1,line=1)

        #text describing the parameters used to generate the plot
        mtext(text=parameters,side=1,line=6,cex=0.7,adj=0.95)
        
        #Close all the opened split.screens and release the figure
        close.screen(all.screens=T)
        dev.off()
    }
}


