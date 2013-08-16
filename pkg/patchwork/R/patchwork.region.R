#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#patchwork.region() ; a function for zooming in on a specific region and showing the regions genes.
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


patchwork.region <- function(CNfile=NULL,chr,region,hg18=F)
{
    if (is.null(CNfile))
    {
        CNfile = list.files(pattern="*_copynumbers.Rdata")
        if (length(CNfile)==0)
        {
        #cat("Could not find <sample>_copynumbers.Rdata in current working directory. Try again and supply patchwork.regions with the CNfile parameter.")
        stop("Could not find <sample>_copynumbers.Rdata in current working directory. Try again and supply patchwork.regions with the CNfile parameter.")
        }
        if(length(CNfile)>1)
        {
        #cat("There seems to be more than one <sample>_copynumbers.Rdata file in the current working directory. Try again and supply patchwork.regions with the CNfile parameter.")
        stop("There seems to be more than one <sample>_copynumbers.Rdata file in the current working directory. Try again and supply patchwork.regions with the CNfile parameter.")
        }
    }
    
    #Load <sample>_copynumbers.Rdata giving us segs,alf and kbsegs
    load(CNfile)

    # #Set working directory and make that directorys name the samplename. If it is too long, shorten it.
    # setwd(directory)
    # subs <- getSubdirs()
    # if (is.null(subs)) 
    # {                                         ## check samples = subdirectories or a single sample = current directory
    #     subs=thisSubdir()
    #     setwd('..')
    # }
    
    #store chr value in nchr because i made the newbie mistake of having a second variable also
    #named chr.
    nchr = chr  
    
    #get sample name
    tmp_name <- strsplit(CNfile,split="_")
    name=tmp_name[[1]][1]
    
    #if sample name is too long, shorten it to 12 characters
    if(nchar(name)>12)
    {
        name = substring(name,1,12)
    }
    
    #If nchr has been given as character
    if(is.character(nchr) == T)
    {
        #If chr is longer than 2 characters ("chr1" but not "12")
        #take out the number and convert it to numeric
        if(nchar(nchr)>2)
        {
            c = strsplit(nchr,"chr")
            c = as.numeric(c[[1]][2])
        }
        #chr has been given as character but without "chr" infront. Convert to numeric.  
        else
        {
            c = as.numeric(nchr) 
        }
    }
    #else chr has been given as a numeric and can be used directly.
    else
    {
        c = nchr
    }
    
    #Rename parameters of <name>_copynumbers so we can re-use other plots as template.
    chr = as.character(segs$chr)
    start = segs$start  
    end = segs$end
    int = segs$median
    ai = segs$ai
    mchr = as.character(kbsegs$chr)
    mpos = kbsegs$pos
    mval = kbsegs$ratio
    schr = as.character(alf$achr) 
    spos = alf$apos
    sval = (1-alf$amin/alf$amax)
    xlim=c(0,2.5)
    ylim=c(0,1)
    
    #check if hg18 is true
    #then load hg18 knownGene list
    if(hg18==T)
    {
        kg = knownGene_hg18
    }
    else
    {
        kg = knownGene
    }
    
    
    #
    #Here is where pretty normal plotting procedure starts
    #
    
    #Get ideogram  #In sysdata.rda
    #ideogram=getIdeogram()
    
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
    pos <- (start+end)/2
    length=end-start

    size=rep(0.5,length(chr))
    size[length>2000000]=1
    size[length>5000000]=1.5
    size[length>10000000]=2
    
    #Select the chromosome
    this <- ideogram[ideogram$c==c,]
    
    #Extract regions start an stop
    #region inputted in form start:stop
    Rstart = min(region)
    Rend = max(region)
    
    #Initialize jpeg
    jpeg(paste(name,'_',this$chr,'_region_',Rstart,"-",Rend,'.jpg',sep=''),width=11.7,height=8.3,units="in",res=300)
    
    #split plot into desired formation
    split.screen(as.matrix(data.frame(left=c(rep(0.05,4),rep(0.47,3)),
                                      right=c(0.45,rep(1,6)),
                                      bottom=c(0.50,0.05,0.2,0.3,0.50,0.70,0.75),
                                      top = c(0.95,0.2,0.3,0.45,0.70,0.75,0.95)))) #screen 1-5
    
    
    #Screen configuration overview
    #-------------------------------------------------------------
    #-------------------------------------------------------------
    #                             |                              |
    #                             |              7               |
    #                             |                              |
    #                             |------------------------------| 
    #            1                |              6               |
    #                             |------------------------------|
    #                             |                              |
    #                             |              5               |
    #                             |                              | 
    #-----------------------------|------------------------------|
    #                                                            |
    #                             4                              |
    #------------------------------------------------------------|
    #                             3                              |
    #------------------------------------------------------------|
    #                                                            |
    #                             2                              |
    #-------------------------------------------------------------
    
    
    #------------------------------------------------------------
    #Left side of the plot (Whole genome)
    #------------------------------------------------------------
    
    #index of overlapping the selected chromosome region to the chr object
    wix <- chr==as.character(this$chr) & (((Rstart <= start) & (Rend >= end)) | ((Rstart >= start) & (Rstart <= end)) | ((Rend <= end) & (Rend >= start)))
    ix <- chr==as.character(this$chr)
    
    #Create an index of colors relating to the positions and lengths on this chromosome
    col <- rep('#B0B0B030',length(chr))       
    col[wix & (pos < this$mid)] <- paste(colors_p(sum(wix & (pos < this$mid))), '70', sep='')  
    col[wix & (pos > this$mid)] <- paste(colors_q(sum(wix & (pos > this$mid))), '70', sep='')
    
    #Go to screen 1
    screen(1)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0.25,0))
    
    #Whole genome overview plot
    #Note that we are plotting chrY in the background, despite the fact that it is not an "active" chromosome.
    plot(c(int[!wix],int[wix]),c(ai[!wix],ai[wix]),
         pch=16,
         cex=c(size[!wix],size[wix]),
         main = "",
         xlab = "",
         ylab = "",
         axes=F, #Remove axis
         col = c(col[!wix],col[wix]),
         xlim = xlim,
         ylim = ylim)
    
    #Insert Y and X axis
    axis(side=2,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),las=1)
    axis(side=1,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=0,to=2.5,by=0.1))
    
    #Titles,date/time and axis labels
    mtext(text="Normalized coverage",side=1,line=1.1,cex=1)
    mtext(text="Allelic imbalance",side=2,line=1.5,cex=1)
    mtext(paste("Detailed view of sample: ",name,"\n Chromosome ",c,", Region: ",Rstart,"-",Rend,sep=""),side=3)
    
    #------------------------------------------------------------
    #Top Right - Signal
    #------------------------------------------------------------
    
    #Go to screen 7
    screen(7)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0,0))
    par(xpd=T)
    
    #Select the correct chromosome and remove stuff lower than -1 
    mix <- mchr==as.character(this$chr) #& mval>(-1)
    
    #Create an index of colors relating to the positions and lengths on this chromosome
    col=rep('#000000',sum(ix))
    col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
    col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))
    
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
    ymin=0
    ymax=2.5
    seqminmax=seq(ymin,ymax,by=0.5)
    
    #Plot coverage over position
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
    
    #Add legend for red region
    mtext(text=expression(bold("Selected region")),side=3,col='#E8000070',cex=0.8,adj=0.05)
    
    #Add colored segments based on the log-ratio data
    segments(x0=start[ix],x1=end[ix],
             y0=int[ix],y1=int[ix],                
             col=col,
             lwd=4)
    
    #Add a rectangle showing which region has been chosen
    rect(xleft=Rstart,xright=Rend,ybottom=ymin,ytop=ymax,col='#E8000040',lty=0)
    
    #Add X and Y axis. The (side=4) axis is just a black line showing where the data ends
    axis(side=2,tck=0.926,col.ticks='#808080',
         #Set axis
         #at=seq(from=-1,to=1.5,by=0.5),
         #labels=c("-1","-.5","0",".5","1","1.5"),
         #Dynamic axis
         at=seqminmax,
         #labels=paste(seqminmax,'',sep=''),
         cex.axis=0.6,pos=0,hadj=1.3,las=1)
    
    axis(side=4,labels=F,tck=0,pos=this$length,
         at=seqminmax)
    
    axis(side=1,tck=0.926,col.ticks='#808080',at=seq(5e6,this$length,by=5e6),
         labels=FALSE,cex.axis=0.6,pos=ymin)
    
    #Add Y axis label & date/time
    mtext("Norm. coverage",side=2,line=0.3)
    mtext(format(Sys.time(),"%Y-%m-%d %H:%M"),side=3,cex=0.6,adj=0.95)
    
    #------------------------------------------------------------
    #Middle Right - Cytobands
    #------------------------------------------------------------
    
    #Go to screen 6
    screen(6)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0,0))
    
    
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
    
    #Add a rectangle showing which region has been chosen
    rect(xleft=Rstart,xright=Rend,ybottom=0,ytop=0.3,col='#E8000040',lty=0)
    
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
    #Bottom Right - Allelic imbalance
    #------------------------------------------------------------ 
    
    #Go to screen 5
    screen(5)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0,0))
    
    
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
    axis(side=2,tck=0.926,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=0,hadj=1.3,las=1)
    axis(side=4,labels=F,tck=0,pos=this$length)
    axis(side=1,tck=0.925,col.ticks='#808080',at=seq(5e6,this$length,by=5e6),
         labels=paste(seq(5,round(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
    
    #Add a rectangle showing which region has been chosen
    rect(xleft=Rstart,xright=Rend,ybottom=0,ytop=1,col='#E8000040',lty=0)
    
    #Add X and Y label
    mtext("Allelic imbalance",side=2,line=0.3)
    mtext("Position (Mb)",side=1,line=1)
    
    #------------------------------------------------------------
    #Top Bottom - Detailed region coverage view
    #------------------------------------------------------------ 
    
    #Go to screen 4
    screen(4)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0,0))
    
    
    #Select the correct chromosome and region and remove stuff lower than -1 
    mix <- mchr==as.character(this$chr) & (mpos >= Rstart) & (mpos <= Rend) #& mval>(-1)
    
    wix <- chr==as.character(this$chr) & (((Rstart <= start) & (Rend >= end)) | ((Rstart >= start) & (Rstart <= end)) | ((Rend <= end) & (Rend >= start)))
    
    cix = ((Rstart >= start) & (Rstart <= end))
    start[cix] = Rstart
    cix = ((Rend <= end) & (Rend >= start))
    end[cix] = Rend
    
    #Predefine ymin ymax and sequence between them
    ymin=0
    ymax=2.5
    seqminmax=seq(ymin,ymax,by=0.5)
    
    #Plot log-ratio over position
    plot(mpos[mix],mval[mix],
         pch=20,
         cex=0.5,
         main = "",
         xlab = "",
         ylab = "",
         axes=F,
         col = '#00000040',
         xlim = c(Rstart,Rend),
         ylim = c(ymin,ymax))#c(round(min(mval[mix]),digits=1),round(max(mval[mix]),digits=1)))#c(-1,1.5))
    
    #Create an index of colors relating to the positions and lengths on this chromosome
    col=rep('#000000',sum(wix))
    col[pos[wix] < this$mid] <- colors_p(sum(pos[wix] < this$mid))
    col[pos[wix] > this$mid] <- colors_q(sum(pos[wix] > this$mid))
    
    #Add colored segments based on the log-ratio data
    segments(x0=start[wix],x1=end[wix],
             y0=int[wix],y1=int[wix],                
             col=col,
             lwd=4)
    
    #Add X and Y axis. The (side=4) axis is just a black line showing where the data ends
    axis(side=2,tck=0.926,col.ticks='#808080',
         #Set axis
         #at=seq(from=-1,to=1.5,by=0.5),
         #labels=c("-1","-.5","0",".5","1","1.5")
         #Dynamic axis
         at=seqminmax,
         #labels=paste(seqminmax,'',sep=''),
         ,cex.axis=0.6,pos=Rstart,hadj=1.3,las=1)
    
    axis(side=4,labels=F,tck=0,pos=Rend)
    #axis(side=1,tck=0.926,col.ticks='#808080',cex.axis=0.6,pos=-1)
    
    #Add X Y axis label
    mtext("Norm. coverage",side=2,line=-0.8)
    #mtext("Position (Mb)",side=1,line=1)
    
    #------------------------------------------------------------
    #Middle Bottom - Region Gene view
    #------------------------------------------------------------ 
    
    #Go to screen 3
    screen(3)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0.25,0))
    
    #par(xpd=T)
    
    #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
    plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,0.48),xlim=c(Rstart,Rend))
    
    kg=kg[order(kg$chr,kg$gtxEnd),]
    
    #Create and index of the genes within region (rstart to rend)
    #gix <- kg$chr==this$chr & (((Rstart <= kg$gtxStart) & (Rend >= kg$gtxEnd)) | ((Rstart >= kg$gtxStart) & (Rstart <= kg$gtxEnd)) | ((Rend <= kg$gtxEnd) & (Rend >= kg$gtxStart)))
    gix <- kg$chr==as.character(this$chr) & ((Rstart <= kg$gtxEnd) & (Rend >= kg$gtxStart)) 
    
    #If they are only partial overlapping, remove the parts outside of the selected region
    cix = ((Rstart >= kg$gtxStart) & (Rstart <= kg$gtxEnd))
    kg$gtxStart[cix] = Rstart
    cix = ((Rend <= kg$gtxEnd) & (Rend >= kg$gtxStart))
    kg$gtxEnd[cix] = Rend
    
    # #Check if any genes are overlapping so they are shown on top of or below eachother
    d1 = c(rep(0.42,length(which(gix))))
    d2 = c(rep(0.36,length(which(gix))))
    j = 1
    
    for(i in min(which(gix)):max(which(gix)))
    {
        #if this gene overlaps with the gene infront of it
        if(kg$gtxEnd[i] >= kg$gtxStart[i+1])
        {
            d1[j+1] = d1[j]-0.06
            d2[j+1] = d2[j]-0.06  
        }
        j = j + 1
    }
    
    #Add the genes to plot as horizontal rectangles that cover the transcribed region
    rect(xleft=kg$gtxStart[gix],xright=kg$gtxEnd[gix],
         #Between 0 - 0.3
         ybottom=d2,ytop=d1,
         col=rgb(kg$R[gix],kg$G[gix],kg$B[gix],maxColorValue=255),lty=1)
    
    #Add legend
    #
    text(x=Rstart+(6.2*(Rend-Rstart)/10),y=0.46,labels=expression(bold("Feature in PDB")),col=rgb(0,0,0,maxColorValue=255),cex=0.6)
    text(x=Rstart+(7.6*(Rend-Rstart)/10),y=0.46,labels=expression(bold("RefSeq,Swissprot or CCDS validated")),col=rgb(12,12,120,maxColorValue=255),cex=0.6)
    text(x=Rstart+(9*(Rend-Rstart)/10),y=0.46,labels=expression(bold("Other RefSeq")),col=rgb(80,80,160,maxColorValue=255),cex=0.6)
    text(x=Rstart+(9.8*(Rend-Rstart)/10),y=0.46,labels=expression(bold("non-RefSeq")),col=rgb(130,130,210,maxColorValue=255),cex=0.6)
    
    #Add Y label
    mtext("Known genes",side=2,line=-1.8,,cex=0.9,las=1)
    
    #legend(c(Rstart,-1),legend=c("Black","Dark blue","Medium blue","Light blue"),cex=0.4,xpd=T,
    #         col=c(rgb(130,130,210,maxColorValue=255),
    #               rgb(12,12,120,maxColorValue=255),
    #               rgb(0,0,0,maxColorValue=255),
    #               rgb(80,80,160,maxColorValue=255)))
    
    for(i in 1:length(kg$gAlias[gix]))
    {
        text(x=(kg$gtxStart[gix][i]+((kg$gtxEnd[gix][i]-kg$gtxStart[gix][i])/2)),
             #This places the gAlias text at same hight level as the transcribed region rectangle.
             #y=((d1[i]+d2[i])/2),
             #This places the gAlias text at the bottom
             y=0.08,
             labels=kg$gAlias[gix][i],
             #This makes the text vertical
             srt=90,
             cex = 0.4,xpd=T)
    }
    
    
    
    #------------------------------------------------------------
    #Bottom Bottom - Detailed region allele frequency view
    #------------------------------------------------------------ 
    
    #Go to screen 2
    screen(2)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0.25,0))
    
    
    #Index of correct chromsome and allele frequency not in either 0 or 1
    six <- schr==as.character(this$chr) & !(sval %in% c(0,1)) & (spos >= Rstart) & (spos <= Rend)
    
    #plot allele frequency over position
    plot(spos[six],sval[six],
         pch=20,
         cex=0.5,
         main = "",
         xlab = "",
         ylab = "",
         #yaxt="n",
         axes=F,
         col = '#00000040',
         xlim = c(Rstart,Rend),
         ylim = c(0,1))
    
    #Add X and Y axis as well as a line to the rightmost of the data (side=4)
    axis(side=2,tck=0.926,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=Rstart,las=1)
    axis(side=4,labels=F,tck=0,pos=Rend)
    #axis(side=1,tck=0.925,col.ticks='#808080',at=seq(5e6,this$length,by=5e6),
    #      labels=paste(seq(5,round(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
    axis(side=1,cex.axis=0.6)
    
    #Add X and Y label
    mtext("Allelic imbalance",side=2,line=-0.8)
    mtext("Position (Mb)",side=1,line=1)
    
    
    #Close all the opened split.screens and release the figure
    close.screen(all.screens=T)
    dev.off()
}






