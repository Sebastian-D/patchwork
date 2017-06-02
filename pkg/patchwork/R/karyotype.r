# karyotype <- function(chr,start,end,int,ai,
# 						name='',xlim=c(0,2.4),ylim=c(0,1))
# 	{
#     png(paste(name,'.karyotype.png',sep=''),width=1300,height=1300)

#     	#packagepath = system.file(package="patchwork")

#     	#load(paste(packagepath,"/data/ideogram.RData",sep=""))

#     	data(ideogram,package="patchworkData")

#     colors_p <- colorRampPalette(c("#2754e8","#ad09e8"),space="rgb")
#     colors_q <- colorRampPalette(c("#c309e8","#e81e2c"),space="rgb")

#     	layout(matrix(1:25,nrow=5,byrow=T), widths=1,heights=1)

#     	aix=ai!=0
#     	chr=chr[aix]
#     	start=start[aix]
#     	end=end[aix]
#     	int=int[aix]
#     	ai=ai[aix]
#     	pos <- (start+end)/2
#     	length=end-start

#     	size=rep(0.5,length(chr))
#     	size[length>2000000]=1
#     	size[length>5000000]=1.5
#     	size[length>10000000]=2


#     	for (c in 1:24)
#     		{

#         	this <- ideogram[ideogram$c==c,]
#         	ix <- chr==this$chr
#         	ix[is.na(ix)]=FALSE

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
#         		pch=16,
#         		cex=c(size[!ix],size[ix]),
#         		cex.lab=3,
#         		mar=c(0.1,0.1,0.1,0.1),
# 				main = "",
# 				xlab = this$chr,
# 				ylab = "",
# 				col = c(col[!ix],col[ix]),
# 				xlim = xlim,
# 				ylim = ylim,
# 		    	yaxt="n",
#         		xaxt="n")
#         	axis(1,cex.axis=1.5)
#         	axis(2,cex.axis=1.5)
#             if(c==1)
#                 {
#                 title(main=paste("Sample: ",name,sep=""),cex.main=2)
#                 }
# 			par(new=F)
#     		}
#     dev.off()
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

#source("/Users/S_D/Documents/svn/patchwork/pkg/patchwork/R/karyotype.r")

#OverviewPlot(as.character(segs$chr),segs$start,segs$end,segs$median,segs$ai,ideogram=ideogram,as.character(kbsegs$chr),kbsegs$pos,kbsegs$ratio,alf$achr,alf$apos,(1-alf$amin/alf$amax),name="wip",chroms=chroms)


karyotype <- function(chr,start,end,int,ai,mchr,mpos,mval,schr,spos,sval,name='',xlim=c(0,2.5),ylim=c(0,1),
    MAPD=" ",MHOF=" ")
{

    #Import ideogram information. Object "chroms" should come from sysdata.rda
    #data(ideogram,package="patchworkData")

    colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")

    ai[is.na(ai)]=0

    aix=ai!=0
    chr=chr[aix]
    start=start[aix]
    end=end[aix]
    int=int[aix]
    ai=ai[aix]
    pos <- (start+end)/2
    length=end-start

    size=rep(0.1,length(chr)); size[int<xlim[1]]=0; size[int>xlim[2]]=0
    #Slightly larger circles for larger segments for (hopefully) clearer clusters
    size[length>2000000]=0.3
    pch=rep(20,length(chr)); pch[chr=='chrX']=4
    #size[length>2000000]=0.4
    #size[length>5000000]=0.6
    #size[length>10000000]=0.8

    #Define the name,dimensions and resolution of the plot.
    jpeg(paste(name,'_overview.jpg',sep=''),width=11.7,height=8.3,units="in",res=300)

    #split the plot into desired formation
    split.screen(figs=c(2,1)) #two rows, one column
    #This details the margins (left,right,bottom,top) of all the screens.
    #See split.screen section of http://seananderson.ca/courses/11-multipanel/multipanel.pdf
    #for the most complete explanation.
    split.screen(as.matrix(data.frame(left=c(rep(seq(from=0.05,to=7*(0.9/8)+0.05,by=(0.9/8)),3)),
                                      right=c(rep(seq(from=(0.9/8)+0.05,to=8*(0.9/8)+0.05,by=(0.9/8)),3)),
                                      bottom=c(rep(2*(0.75/3)+0.15,8),rep((0.75/3)+0.15,8),rep(0.15,8)),
                                      top=c(rep(0.9,8),rep(2*(0.75/3)+0.15,8),rep((0.75/3)+0.15,8)))),1) #screen 1 becomes screen 3 - 26

    split.screen(as.matrix(data.frame(left=c(rep(0.05,3)),
                                      right=c(rep(1,3)),
                                      bottom=c(0.15,0.5,0.55),
                                      top=c(0.5,0.55,1))),2) #screen 2 become screen 27-29

    #Screen configuration overview
    #------------------------------------------------------------
    #|   3  |   4    |   5   |  6   |   7  |  8   |  9   |  10  |
    #------------------------------------------------------------
    #|  11  |   12   |   13  |  14  |   15 |  16  |  17  |  18  |
    #------------------------------------------------------------
    #|  19  |   20   |   21  |  22  |   23 |  24  |  25  |  26  |
    #------------------------------------------------------------
    #|                            29                            |
    #------------------------------------------------------------
    #|                            28                            |
    #------------------------------------------------------------
    #|                            27                            |
    #------------------------------------------------------------


    #------------------------------------------------------------
    #Top part of the plot
    #------------------------------------------------------------

    #nchr=24; if (sum(chr=='chrY')==0) nchr=23
    #Loop over and plot the 24 chromosomes
    for (c in 1:24)
    {
        #Pick a chromosome
        this <- ideogram[ideogram$c==c,]
        #Extract that chromosomes information
        ix <- chr==as.character(this$chr)
        #Add color depending on the length of this chromosome
        col <- rep('#B0B0B030',length(chr))
        col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')
        col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='')

        #Enter different screens to add whole genome plots
        screen(c+2)

        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        par(xpd=NA)

        #Plot whole genome with chromosome in question colored
        plot(c(int[!ix],int[ix]),c(ai[!ix],ai[ix]),
             pch=c(pch[!ix],pch[ix]),
             cex=c(size[!ix],size[ix]),
             main = "",
             xlab = '',
             ylab = "",
             col = c(col[!ix],col[ix]),
             xlim=xlim,ylim=c(0,1),
             axes=F)

        #grid lines
        segments(
            x0=c(0.5,1,1.5),x1=c(0.5,1,1.5),
            y0=0,y1=1,
            col='#D3D3D360',
            lwd=1)

        #Add boxes around the plot to distinguish the chromosomes further
        box()

        #Add chromosome number inside plot
        if(c == 23)
        {
            mtext("X", side = 3, line = -1, adj = 0.8, cex = 1)
        } else if(c == 24)
        {
            mtext("Y", side = 3, line = -1, adj = 0.8, cex = 1)
        } else
        {
            mtext(c, side = 3, line = -1, adj = 0.8, cex = 1)
        }

        #If we are one of the chromosomes plotted to the left edge or bottom, add axis
        if(c %in% c(1,9,17)) axis(side=2,cex.axis=0.6,tck=-0.05,at=seq(from=0,to=0.8,by=0.2),las=1)
        if(c %in% c(seq(from=17,to=24,by=1))) axis(side=1,cex.axis=0.6,at=seq(from=0,to=2,by=0.5))
        if(c %in% c(8,16,24)) axis(side=4,cex.axis=0.6,tck=-0.05,at=seq(from=0,to=0.8,by=0.2),las=1)

        #Add titles with sample name (name of folder containing sample) and date
        if(c==4)
        {
            mtext(name,side=3)
        }
        if(c==8)
        {
            mtext(paste('MAPD: ',MAPD,', MHOF: ',MHOF,'\n',format(Sys.time(),"%Y-%m-%d %H:%M")),side=3,cex=0.6,adj=0.95)
        }
        if(c==9)
        {
            mtext("Allelic imbalance",side=2,line=1.5)
        }
        if(c==20)
        {
            mtext("Normalized coverage",side=1,line=1.5)
        }
    }

    #------------------------------------------------------------
    #Middle - Signal
    #------------------------------------------------------------

    #Go to screen 29 for plotting
    screen(29)

    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))

    #Calculate previous distance of whole genome so the next chromosome is added
    #at the correct coordinate
    pre=rep(NA,24)
    for (c in 1:24)
    {
        this <- ideogram[ideogram$c==c,]
        ix <- chr==as.character(this$chr)
        mix <- mchr==as.character(this$chr) #& mval>=(-2)
        six <- schr==as.character(this$chr)

        #Predefine ymin ymax and sequence between them
        # ymin=floor(2*min(int,na.rm=T))/2
        # if(ymin > -1)
        # {
        #     ymin = -1
        # }
        # ymax=ceiling(2*max(int,na.rm=T)+1)/2
        # if(ymax < 1)
        # {
        #     ymax = 1
        # }
        ymin=0
        ymax=2.5
        seqminmax=seq(ymin,ymax,by=0.5)

        if(c>1)
        {
            prelength=prelength+ideogram$length[ideogram$c==(c-1)]
        }
        else
        {
            prelength=0
        }

        #position used for signal.
        mpos[mix] = mpos[mix] + prelength
        #Start and end used for colored segment placement (Adding whole preceding genome to each position)
        start[ix] = start[ix] + prelength
        end[ix] = end[ix] + prelength
        #pre used to position chromosome markers. Seen bottom of plot.
        pre[c] = prelength+(this$length/2)
        #spos used later for allelic imbalance plot. Same principle as mpos.
        spos[six] = spos[six] + prelength
    }

    #Only plot mval higher than -2
    #ix=mval>=(-2)
    #mpos = mpos[ix]
    #mval = mval[ix]

    #Remove every third value to make plotting more sparse in this overview. (saves space and time)
    ix = seq(3,to=length(mpos),by=3)
    mpos = mpos[-ix]
    mval = mval [-ix]

    #Plot signal over complete genome
    plot(mpos,mval,
         pch=20,
         cex=0.25,
         main='',
         xlab = "",
         ylab = "",
         col = '#00000010',
         xaxt="n",
         axes=F,
         ylim = c(ymin,ymax),
         xlim = c(0,max(mpos))
    )

    #Add axis to the left & right of signal
    axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
    axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=max(mpos),las=1)
    mtext("Normalized coverage",side=2,line=0)

    #Add grey segments
    segments(
        y0=c(0.5,1,1.5),y1=c(0.5,1,1.5),
        x0=0,x1=max(mpos),
        col='#D3D3D360',
        lwd=1)

    #Add colored segment information for each chromosome, red to blue gradient.
    for(c in 1:24)
    {
        this <- ideogram[ideogram$c==c,]
        ix <- chr==as.character(this$chr)
        col <- rep('#B0B0B030',length(chr))
        col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')
        col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='')
        col=rep('#000000',sum(ix))
        col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
        col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))

        #Add segments
        segments(x0=start[ix],x1=end[ix],
                 y0=int[ix],y1=int[ix],
                 col=col,
                 lwd=2)
    }

    #Add a bar between chromosomes to distinguish them
    segments(
        x0=c(chroms$before,sum(chroms$length),max(mpos)),x1=c(chroms$before,sum(chroms$length),max(mpos)),
        y0=-100,y1=100,
        col='#000000',
        lwd=1)
    #Top and bottom bar
    #segments(
    #    x0=c(0,0),x1=c(max(mpos),max(mpos)),
    #    y0=c(ymin,ymax),y1=c(ymin,ymax),
    #    col='#000000',
    #    lwd=1)

    #------------------------------------------------------------
    #Cytobands
    #------------------------------------------------------------

    #Go to screen 28 for plotting
    screen(28)

    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))

    #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
    plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,0.3),xlim=c(0,max(mpos)))


    mtext(text='',side=2,las=1,line=-1.2)

    #Add cytoband information as differently colored rectangles
    rect(xleft=chromData$genomeStart,xright=chromData$genomeEnd,
         ybottom=0.15-chromData$thickness*0.03,ytop=0.15+chromData$thickness*0.03,
         col=paste(chromData$col,'80',sep=''),
         lty=0)

    #Add a bar between chromosomes to distinguish them
    segments(
        x0=c(chroms$before[-1],sum(chroms$length)),x1=c(chroms$before[-1],sum(chroms$length)),
        y0=-100,y1=100,
        col='#000000',
        lwd=1)

    #------------------------------------------------------------
    #Bottom - AI
    #------------------------------------------------------------

    #Go to screen 27
    screen(27)

    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))

    #Index to avoid spos/sval values that are all the way at 0 or 1.
    ix = !(sval %in% c(0,1))

    #Plot allelic imbalance over whole genome
    plot(spos[ix],sval[ix],
         pch=20,
         cex=0.35,
         cex.axis=1,
         main='',
         xlab = "",
         ylab = "",
         col = '#00000003',
         xaxt="n",
         axes=F,
         ylim = c(0,1),
         xlim = c(0,max(mpos))
    )

    #Add axis to the left,right and below of AI. The below axis is the chromosome numbers 1-24.
    axis(side=2,tck=-0.04,at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=0,las=1)
    axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X","Y"),cex.axis=0.55,lty=0)#,tck=0,col.ticks='#00000000')
    axis(side=4,tck=-0.04,at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=max(mpos),las=1) #
    mtext("Allelic imbalance",side=2,line=0)
    mtext("Chromosomes",side=1,line=1.5,adj=0.4)

    #Add a bar between chromosomes to distinguish them
    segments(
        x0=c(chroms$before,sum(chroms$length),max(mpos)),x1=c(chroms$before,sum(chroms$length),max(mpos)),
        y0=0,y1=100,
        col='#000000',
        lwd=1)


    #Close all the opened split.screens and release the figure
    close.screen(all.screens=T)
    dev.off()

    file.copy(paste(name,'_overview.jpg',sep=''),paste("../",name,'_overview.jpg',sep=''),overwrite=T)

}


#test
