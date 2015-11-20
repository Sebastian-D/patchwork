## load('~/patchwork/pkg/TAPS/R/sysdata.rda')
## load('shortRegions.Rdata')
## load('allRegions.Rdata')

## TAPS source code reimplemented 2013 Markus (Rasmussen) Mayrhofer, Sebastian DiLorenzo

#########    ###    ########   ######     ########  ##        #######  ######## 
#   ##      ## ##   ##     ## ##    ##    ##     ## ##       ##     ##    ##    
#   ##     ##   ##  ##     ## ##          ##     ## ##       ##     ##    ##    
#   ##    ##     ## ########   ######     ########  ##       ##     ##    ##    
#   ##    ######### ##              ##    ##        ##       ##     ##    ##    
#   ##    ##     ## ##        ##    ##    ##        ##       ##     ##    ##    
#   ##    ##     ## ##         ######     ##        ########  #######     ##  

TAPS_plot <- function(#samples='all',
                     directory=NULL,autoEstimate=FALSE,
                      bin=250,cores=1,matched=FALSE,allelePeaks=FALSE) {
    # cores=1
    #Automatically check, and if needed install, packages stats and fields
    
    #Load stats. It should be in all, at least semi-new, R distributions so we dont need to install.package it or
    #pre-install it
    #library(stats)
    suppressPackageStartupMessages(library(stats))
    suppressPackageStartupMessages(library(DNAcopy))
    suppressPackageStartupMessages(library(fields))
    suppressPackageStartupMessages(library(xlsx))    
    suppressPackageStartupMessages(library(foreach))
    try(suppressPackageStartupMessages(library(doMC)),silent=T)
    try(suppressPackageStartupMessages(registerDoMC(cores=cores)),silent=T)


    #list.of.packages <- c("stats", "fields")
    #new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    #if(length(new.packages)) install.packages(new.packages)
    
    if (is.null(directory))
    {
        print("Using working directory")
        directory = "."
        #cat("You have not assigned a directory containing one or more folders of samples for TAPS_plot to execute. \n")
        #cat("Example: \"/user/mysamples/\" or, to run it in your current working directory, \".\" \n")
        #directory = readline("Please supply such a directory now: ")
    }
    
    ## This function takes a directory as input, then builds short-segment TAPS scatter plots for each sample (subdirectory) in the directory.
    setwd(directory)
    subs <- getSubdirs()
    subs=subs[subs!='frequencies' & subs!='frequencies_comp']
    subsToSampleData=NULL
    if (is.null(subs)) {                                         ## check samples = subdirectories or a single sample = current directory
        subs=thisSubdir()
        setwd('..')
        subsToSampleData <- getSubdirs()
    }
    
   # create SampleData file if there is none.   
    if (length(grep('SampleData.csv',dir()))==0) {
        if(!is.null(subsToSampleData))
            {
            sampleData <- data.frame(Sample=subsToSampleData,cn1= NA, cn2=NA, cn3=NA, loh=NA, MAPD=NA, MHOF=NA)
            }
        else
            {
            sampleData <- data.frame(Sample=subs,cn1= NA, cn2=NA, cn3=NA, loh=NA, MAPD=NA, MHOF=NA)
            }
        save.txt(sampleData,'SampleData.csv')
    } else {
        sampleData=load.txt('SampleData.csv')
    }

    #if (samples[1]=='all') samples=rep(T,length(subs))
    #if (is.logical(samples)) samples=which(samples)
    #subs=subs[samples]
    root <- getwd()
    print(paste('root: ',root,sep=''))
    
    # for (i in 1:length(subs)) {
    #junk only stores the list from foreach.
    sampleData$MHOF <- sampleData$MAPD <- NA
    junk <- foreach (i = 1:length(subs)) %dopar% {
        setwd(paste(root,'/',subs[i],sep=''))
        name <- subs[i]
        #if (length(grep('sampleData.csv',dir()))==0) save.txt(data.frame(Sample=name,cn2='',delta='',loh='',completed.analysis='no'),file='sampleData.csv')
        # cat(' ..loading', subs[i])
        # print(paste(i,': ',subs[i],': Opening',sep=''))
        print(paste(i,'/',length(subs),': ',subs[i],' Loading',sep=''))
        
        if(length(grep("*.cychp",dir()))==1)				##cychp sample
        {
            ##read CYCHP file from ChAS with log-ratio info##
            ################################################
            suppressPackageStartupMessages(library(affxparser))
            
            name=list.files(".",pattern="*.cychp")
            temp=readCcg(name)
            tempLog2=temp$dataGroups$ProbeSets$dataSets$CopyNumber$table[,1:4]
            #Value=NULL
            Start=tempLog2$Position
            End=tempLog2$Position
            Value=tempLog2$Log2Ratio
            tempChr=tempLog2$Chromosome
            nrows=dim(tempLog2)[1]
            Chromosome=rep(NA,nrows)
            for (i in 1:nrows) {
                Chromosome[i]=paste('chr', tempChr[i], sep="")
            }
            
            Log2=as.data.frame(cbind(Chromosome, Start, End, Value))
            
            colnames(Log2)=c("Chromosome","Start","End","Value")
            
            Log2$Chromosome = as.character(Log2$Chromosome)
            Log2$Chromosome[Log2$Chromosome=="chr24"] = "chrX"
            Log2$Chromosome[Log2$Chromosome=="chr25"] = "chrY"
            Log2$Chromosome = as.factor(Log2$Chromosome)
            
            Log2$Start = as.integer(as.character(Log2$Start))
            Log2$End = as.integer(as.character(Log2$End))
            Log2$Value = as.numeric(as.character(Log2$Value))
            
            ##read txt file from ChAS with snp info##
            #########################################
            
            name=c(list.files(".",pattern="*.cy750K.txt"),list.files(".",pattern="*.cyhd.txt"))
            tempalf=read.table(file=name, header=FALSE, sep="\t", skip=12)
            SignalA=tempalf$V4
            SignalB=tempalf$V5
            tempChr=tempalf$V8
            Start=tempalf$V9
            End=tempalf$V9
            Chromosome=Value=rep(NA,nrows)
            nrows=dim(tempalf)[1]
            for (i in 1:nrows) {
                Value[i]=SignalB[i]/(SignalB[i]+SignalA[i])
                Chromosome[i]=paste('chr', tempChr[i], sep="")
            }
            
            alf=as.data.frame(cbind(Chromosome, Start, End, Value))
            
            colnames(alf)=c("Chromosome","Start","End","Value")
            
            alf$Start = as.integer(as.character(alf$Start))
            alf$End = as.integer(as.character(alf$End))
            alf$Value = as.numeric(as.character(alf$Value))
        } else if (file.exists('rawcopy.Rdata')) {
            load('rawcopy.Rdata')
            Log2=probes.txt[,2:5]
            alf=snps.txt[,2:5]
        } else {
            Log2 <- readLog2()                                     ## Log-R
            alf <- readAlf()                                         ## Allele Frequency
        }
        
        #browser()
        
        Log2=Log2[!is.nan(Log2$Value),]
        Log2=Log2[!is.na(Log2$Value),]
        Log2 <- Log2[which(Log2$Value != -Inf & Log2$Value != +Inf ),]
        Log2=chromFix(Log2)
        alf=alf[!is.nan(alf$Value),]
        alf=alf[!is.na(alf$Value),]
        alf$Value[alf$Value<0]=0; alf$Value[alf$Value>1]=1
        alf <- alf[which(alf$Value != -Inf & alf$Value != +Inf ),]
        alf=chromFix(alf)
        
        segments <- readSegments()                                ## segments if available (CBS recommended)
 
        # cat(' ..processing')
        if (is.null(segments)) {                                 ## segmentation using DNA-copy if needed (must then be installed)
            segments <- segment_DNAcopy(Log2)
            save.txt(segments,'_segments.txt') 
        }

        #Remove NA values that some samples give.
        segments <- segments[!is.nan(segments$Value),]
        segments <- segments[!is.na(segments$Value),]
        
        segments$Value <- segments$Value-median(Log2$Value,na.rm=T)     ## Median-centering
        Log2$Value <- Log2$Value-median(Log2$Value,na.rm=T)             ## Median-centering
        
        #browser()
        allRegions=NULL; #if ('allRegions.Rdata' %in% dir()) load('allRegions.Rdata')
        if (is.null(allRegions)) allRegions <- makeRegions(Log2, alf, segments,matched=matched,min=30,allelePeaks=allelePeaks)            ## Calculates necessary data for segments (all functions are in this file)
        save(allRegions,file='allRegions.Rdata')
        regs=NULL;# if ('shortRegions.Rdata' %in% dir()) load('shortRegions.Rdata')
        #browser()
        if (is.null(regs)) {
            regs <- regsFromSegs(Log2,alf,segments,bin=bin,min=30,matched=matched,allelePeaks=allelePeaks)    ## Calculates the same data for shortened segments
            save(regs,file='shortRegions.Rdata')
        }
        
        ## Sample QC 
        # print(paste(i,': ',segments,', ',' sample QC',sep=''))
        # sampleData$MAPD[i] <- MAPD <- round(median(abs(diff(Log2$Value[Log2$Chromosome=='chr1'][order(Log2$Start[Log2$Chromosome=='chr1'])]))),2)
        # sampleData$MHOF[i] <- MHOF <- round(100*median(0.5+abs(0.5-alf$Value[alf$Chromosome!='chrX'])),1)
        # print(round(median(abs(diff(Log2$Value[Log2$Chromosome=='chr1'][order(Log2$Start[Log2$Chromosome=='chr1'])]))),2))
        # print(round(100*median(0.5+abs(0.5-alf$Value[alf$Chromosome!='chrX'])),1))

        sampleData$MAPD[i] <- round(median(abs(diff(Log2$Value[Log2$Chromosome=='chr1'][order(Log2$Start[Log2$Chromosome=='chr1'])]))),2)
        sampleData$MHOF[i] <- round(100*median(0.5+abs(0.5-alf$Value[alf$Chromosome!='chrX'])),1)

            #round(100*median(0.5+regs$hom[regs$scores<0.5],na.rm=T),1)        
        #MAID=round(median(abs(diff(regs$scores[!is.na(regs$scores)]))),3)
        
        #Save for TAPS_region()
        # print(paste(i,': taps_region()',sep=''))
        save(regs,Log2,alf,segments,file="TAPS_plot_output.Rdata")
        
        #save.txt(sampleData,file='../sampleData.csv')
        
        #Clear warnings generated previously so hopefully I can see what is actually causing the program to fail.
        #assign("last.warning", NULL, envir = baseenv())
        
        #Test if hg18 or hg19 should be used. length of (hg18 chr19) > (hg19 chr19)
        hgtest=regs[regs$chr=="chr19",]
        if(hgtest$end[length(hgtest$chr)] > 60000000) {
            hg18=T
        } else {
            hg18=F
        }

        # cat('..plotting.\n')
        OverviewPlot(regs$chr,regs$start,regs$end,regs$logs,regs$scores,hg18=hg18,
                     as.character(Log2$Chromosome),Log2$Start,Log2$Value,as.character(alf$Chromosome),alf$Start,alf$Value,
                     name=name,MAPD=sampleData$MAPD[i],MHOF=sampleData$MHOF[i])                
        
        ## Chromosome-wise plots for manual analysis
        regions=allRegions$regions
        
        karyotype_chroms(regs$chr,regs$start,regs$end,regs$logs,regs$scores,hg18=hg18,
                         as.character(Log2$Chromosome),Log2$Start,Log2$Value,as.character(alf$Chromosome),alf$Start,
                         alf$Value,name=name,MAPD=sampleData$MAPD[i],MHOF=sampleData$MHOF[i])
        
        ## Finally add estimates if needed:
        e=NULL
        if (is.logical(autoEstimate)) {
            if (length(autoEstimate)>1) {
                if (autoEstimate[i]) e=getEstimates(regs$logs,regs$scores)
            } else if (autoEstimate[1]) e=getEstimates(regs$logs,regs$scores)
        } else if (is.numeric(autoEstimate)) { 
            if (i %in% autoEstimate) e=getEstimates(regs$logs,regs$scores)
        } else if (is.character(autoEstimate)) if (name %in% autoEstimate)
            e=getEstimates(regs$logs,regs$scores)
        if (!is.null(e)) {
            sampleData$cn1=e[1]
            sampleData$cn2=e[2]
            sampleData$cn3=e[3]
            sampleData$loh=e[4]
        }
        
        # cat('..done\n')
        # print(paste(i,': ',subs[i],': OK',sep=''))
        print(paste(i,'/',length(subs),': ',subs[i],' OK',sep=''))
        setwd(root)
        1
    }
    save.txt(sampleData,'SampleData.csv')
}
###

#########    ###    ########   ######      ######     ###    ##       ##       
#   ##      ## ##   ##     ## ##    ##    ##    ##   ## ##   ##       ##       
#   ##     ##   ##  ##     ## ##          ##        ##   ##  ##       ##       
#   ##    ##     ## ########   ######     ##       ##     ## ##       ##       
#   ##    ######### ##              ##    ##       ######### ##       ##       
#   ##    ##     ## ##        ##    ##    ##    ## ##     ## ##       ##       
#   ##    ##     ## ##         ######      ######  ##     ## ######## ######## 

###
TAPS_call <- function(samples='all',directory=getwd(),cores=1) {
    # cores=1
    minseg=1
    maxCn=12
    #suppressPackageStartupMessages(library(xlsx))    
    suppressPackageStartupMessages(library(foreach))
    try(suppressPackageStartupMessages(library(doMC)),silent=T)
    try(suppressPackageStartupMessages(registerDoMC(cores=cores)),silent=T)
    
    
    ## TAPS_call outputs the total and minor allele copy numbers of all segments as a text file, and as images for visual confirmation.
    ## sampleInfo_TAPS.txt must be present in each sample folder. If TAPS_plot could not make a good guess of the Log-R of copy number 2 
    ## and the Log-R difference to a deletion, you must interpret the scatter plots and edit sampleInfo_TAPS.txt.
    if (is.null(directory))
    {
        print("No directory supplied, using working directory.")
        directory = "."
        #cat("You have not assigned a directory containing one or more folders of samples for TAPS_call to execute. \n")
        #cat("Example: \"/user/mysamples/\" or, to run it in your current working directory, \".\" \n")
        #directory = readline("Please supply such a directory now: ")
    }
    
    setwd(directory)
    #subs <- getSubdirs()
    
    if (length(grep('SampleData.csv',dir()))==0)
    {
        if (length(grep('SampleData.xlsx',dir()))==1)
        {
            try(library(xlsx))
            sampleData=read.xlsx('SampleData.xlsx',1)
            save.txt(sampleData,'SampleData.csv')
        }
        sampleData=load.txt('SampleData.csv')
    }
    
    if (length(grep('SampleData.csv$',dir()))==1)
    {
        sampleData=load.txt('SampleData.csv')
    } else
    {
        sampleData <- load.txt('../SampleData.csv')
    }
    subs=as.character(sampleData$Sample)
    
    if (is.null(subs)) {
        subs=thisSubdir()
        subs=subs[subs!='frequencies' & subs!='frequencies_comp']
        setwd('..')
    }
    
    if (samples[1]=='all') samples=rep(T,length(subs))
    if (is.logical(samples)) samples=which(samples)
    subs=subs[samples]
    sampleData <- sampleData[samples,]

    # for (i in 1:length(subs)) {
    #junk only stores the list from foreach.
    junk <- foreach(i=1:length(subs)) %dopar% {
        setwd(subs[i])
        name <- subs[i]
        sampleInfo <- sampleData[i,c('cn1','cn2','cn3','loh')]
        if (nrow(sampleInfo)==1) if (sum(is.na(sampleInfo))<4) {
            
                # cat(' ..loading', subs[i])
            print(paste(i,'/',length(subs),': ',subs[i],' Loading',sep=''))
            if (file.exists('rawcopy.Rdata')) {
                load('rawcopy.Rdata')
                Log2=probes.txt[,2:5]
                alf=snps.txt[,2:5]
            } else {
                Log2 <- readLog2()
                alf <- readAlf(localDir)
            }
            segments <- readSegments()
            
            #Some samples contain NA values, we simply remove these.
            Log2=Log2[!is.nan(Log2$Value),]
            Log2=chromFix(Log2)
            
            Log2=Log2[!is.na(Log2$Value),]
            
            alf=alf[!is.nan(alf$Value),]
            alf=alf[!is.na(alf$Value),]
            alf=chromFix(alf)
            
            
            segments <- segments[!is.nan(segments$Value),]
            segments <- segments[!is.na(segments$Value),]    
            segments=chromFix(segments)
            
            
            segments$Value <- segments$Value-mean(Log2$Value) 
            Log2$Value <- Log2$Value-mean(Log2$Value)
            
            # cat(' ..processing.\n')
            
            load('allRegions.Rdata')                            ## These were prepared in TAPS_plot
            load('shortRegions.Rdata')
            #allRegions <- makeRegions(Log2, alf, segments)
            
            ## estimates the Log-R and Allelic Imbalance Ration of all variants up to maxCn
            t <- findCNs(Log2,alf,allRegions,regs,name,maxCn=maxCn,ceiling=1,sampleInfo=sampleInfo) 
            save(t,file='t.Rdata')
            
            u <- setCNs(allRegions,t$int,t$ai,t$model,maxCn)            ## Assigns copy number variant for all segments
            allRegions$regions <- u$regions
            ## adjacent segments with idendical copy number are NOT... merged (except over centromere) and all are saved to a text file
            save.txt(u$regions,file=paste(name,'_segmentCN.txt',sep='')) 
            regions=allRegions$regions

            #save(u$model,file="model.Rdata")
            write.table(t(as.data.frame(u$model)),file='model.txt',row.names=T)

            save(t,regions,file="regions_t.Rdata")


            #Test if hg18 or hg19 should be used. length of (hg18 chr19) > (hg19 chr19)
            hgtest=regions[regions$Chromosome=="chr19",]
            if(hgtest$End[length(hgtest$Chromosome)] > 60000000)
            {
                hg18=T
            } else {
                hg18=F
            }

            #save parameters as strings
            parameters=paste("Parameters given: cn2:",sampleInfo$cn2," delta:",sampleInfo$delta," loh:",sampleInfo$loh)
            
            #karyotype_check(regions$Chromosome,regions$Start,regions$End,regions$log2,regions$imba,regions$Cn,regions$mCn,t,name=name)
            
            karyotype_chromsCN(regions$Chromosome,regions$Start,regions$End,regions$log2,
                               regions$imba,regions$Cn,regions$mCn,hg18=hg18,
                               as.character(Log2$Chromosome),Log2$Start,Log2$Value,as.character(alf$Chromosome),
                               alf$Start,alf$Value,t,name=name,xlim=c(-1,1),ylim=c(0,1),parameters=parameters)
            
            # cat('..done\n')
            print(paste(i,'/',length(subs),': ',name,' ', 'OK',sep=''))
        } else print(paste(i,'/',length(subs),': ',name,' ', 'Skipped',sep=''))
        
        setwd('..')
    }
    #save.txt(sampleData,file='sampleData.csv')
    
}
###
segMatch <- function(segments.txt,probes.txt) {
    ## Returns a data frame with probes.startIx probes.endIx
    #verry fast
    colnames(segments.txt)[colnames(segments.txt)=='chr']='Chromosome'
    colnames(segments.txt)[colnames(segments.txt)=='start']='Start'
    colnames(segments.txt)[colnames(segments.txt)=='end']='End'
    
    segStarts=data.frame(n.seg=1:nrow(segments.txt),n.pos=NA,chr=chrom_n(segments.txt$Chromosome),pos=segments.txt$Start,type='start')
    segEnds=data.frame(n.seg=1:nrow(segments.txt),n.pos=NA,chr=chrom_n(segments.txt$Chromosome),pos=segments.txt$End,type='end')
    lpos=data.frame(n.seg=NA,n.pos=1:nrow(probes.txt),chr=chrom_n(probes.txt$Chromosome),pos=(probes.txt$Start+probes.txt$End)/2,type='logr')
    dummy1=data.frame(n.seg=NA,n.pos=NA,chr=-Inf,pos=-Inf,type='logr')
    dummy2=data.frame(n.seg=NA,n.pos=NA,chr=Inf,pos=Inf,type='logr')
    
    table=rbind(segStarts,segEnds,lpos,dummy1,dummy2)
    table=table[order(table$chr,table$pos),]
    
    start=(table$type=='start')
    end=(table$type=='end')
    #l=(table$type=='logr')
    
    # omitting segment ends, let segment start "n.pos" be that of the next marker
    table$n.pos[!end][start[!end]] <- table$n.pos[!end][which(start[!end])+1]
    
    # omitting segment starts, let segment end "n.pos" be that of the previous marker
    table$n.pos[!start][end[!start]] <- table$n.pos[!start][which(end[!start])-1]
    
    # re-create a data frame with startIx and endIx (in probes.txt) per segment
    indices=data.frame(startIx=table$n.pos[start],endIx=table$n.pos[end])
    
    na=is.na(indices$startIx+indices$endIx)
    na[!na] = indices$startIx[!na]>indices$endIx[!na]
    indices$startIx[na] <- indices$endIx[na] <- NA
    
    return(indices)
}
regsFromSegs <- function (Log2,alf, segments, bin=150,min=1,matched=F,allelePeaks=FALSE) {
    ## This function builds short segments and calcualtes their average Log-R and Allelic Imbalance Ratio.
    rownames(Log2)=1:nrows(Log2)
    rownames(alf)=1:nrows(alf)
    regs=NULL #list('chr'=NULL,'start'=NULL,'end'=NULL,'logs'=NULL,'scores'=NULL,'probes'=NULL,'snps'=NULL)
                                #,'key1'=rep(NA,nrow(Log2)),'key2'=rep(NA,nrow(alf)))
    n=nrow(segments)
    s_check=NULL
    
    pMatch=segMatch(segments,Log2)
    sMatch=segMatch(segments,alf)
    
    for (c in 1:n) { ## for every segment
        #tlog=Log2[Log2$Chromosome==segments$Chromosome[c],] ## Log-R on this chromosome
        #talf=alf[alf$Chromosome==segments$Chromosome[c],] ## Allele Freq on this chromosome
        #tlog=tlog[(tlog$Start>=segments$Start[c])&(tlog$Start<segments$End[c]),] ## Log-R on this segment
        #talf=talf[(talf$Start>=segments$Start[c])&(talf$Start<segments$End[c]),] ## Allele Freq on this segment
        
        tsnps=sMatch[c,2]-sMatch[c,1]+1; if (is.na(tsnps)) tsnps=0   #nrow(talf)    ## number of snps and probes in this segment
        tprobes=pMatch[c,2]-pMatch[c,1]+1; if (is.na(tsnps)) tsnps=0   #nrow(tlog)
        tnregs=max(trunc(tsnps/bin),1) ## Further split into how many (a minimum of 1) subsegments
        tcuts=segments$Start[c] ## The first start pos
        tlength=segments$End[c]-segments$Start[c]    ## Length of this whole segment
        for (i in 1:tnregs) tcuts = c(tcuts, tcuts[1]+round(i/tnregs*tlength)) ## Break the segment at these positions
        
        tab=data.frame(chr=segments$Chromosome[c],start=tcuts[-length(tcuts)],end=tcuts[-1],logs=NA,scores=NA,probes=0,snps=0)
        if (is.null(regs)) regs=tab else regs=rbind(regs,tab)
    }
    
    pMatch=segMatch(regs,Log2)
    sMatch=segMatch(regs,alf)
        
    for (i in 1:nrow(regs)) {    ## build the subsegments
        #regs$chr=c(regs$chr,as.character(segments$Chromosome[c]))    ## Chromosome
        #s_=tcuts[r]                                                     ## Start
        #e_=tcuts[r+1]                                                 ## End
        #thisalf=talf[(talf$Start>=s_)&(talf$Start<=e_),]             ## get the Log-R values
        #thislog=tlog[(tlog$Start>=s_)&(tlog$Start<=e_),]             ## and the allele frequency
        #regs$key1[as.integer(rownames(thislog))]=length(regs$chr)    ## store their positions for fast access during plotting
        #regs$key2[as.integer(rownames(thisalf))]=length(regs$chr)    ## --"--
        if (!is.na(pMatch[i,1])) { 
            regs$logs[i]=median(Log2$Value[pMatch[i,1]:pMatch[i,2]],na.rm=T)  ## store average log ratio of this segment
            regs$probes[i]=pMatch[i,2]-pMatch[i,1]+1                   ## store number of probes
        }
        if (!is.na(sMatch[i,1])) { 
            temp = NA
            try(temp <- allelicImbalance(alf$Value[sMatch[i,1]:sMatch[i,2]],min,matched,allelePeaks),silent=T)
            regs$scores[i]=temp
            regs$snps[i]=sMatch[i,2]-sMatch[i,1]+1                       
        }
    }

    #regs=as.data.frame(regs)
    regs=regs[!is.na(regs$logs),]  ### MODDAT MARKUS MAJ 2013
    return (regs)
}
allelicImbalance <- function (data,min=30,matched=F,allelePeaks=F) {
    #browser()
    if(matched == T ) {
         return(2*median(abs(data-.5),na.rm=T))
    }
    if(allelePeaks == F) {
        if (length(data)>min) {                                    ## Time to calculate Allelic Imbalance Ratio (if enough SNPs)
            t1=sort( abs(data-0.5) )                            ## distance from middle (het) in the allele freq pattern, ascending
            if (length(unique(t1))>3) {                                ## do not attempt clustering with too few snps
                xx=NULL
                try(xx <- kmeans(t1, 2),silent=T)                            ## Attempt k-means (Hartigan-Wong: has proven very stable)
                if (!is.null(xx)) if (min(xx$size) > 0.05*max(xx$size)) {    ## On some occations data quality is poor, requiring 5%+ heterozygous SNPs avoids most such cases.
                    xx=xx$centers
                } else xx=median(t1,na.rm = T)
            } else xx=median(t1,na.rm = T)      
        } else xx=NA
        #try (if (is.na(xx)) xx=0:1, silent=T)
        #try (if (length(xx)==0) xx=0:1, silent=T)
        #regs$scores=c(regs$scores, min(xx)/max(xx) )                ## Allelic Imbalance Ratio = inner / outer cluster.
        #regs$het=c(regs$het, min(xx))                                ## $het and $hom are no longer in use.
        #regs$hom=c(regs$hom, max(xx))
        return( min(xx)/max(xx) )
    }
    #The new shiet
    if(allelePeaks == T) {
        return(1-max(kmeans(data,2)$centers))
    }
}
###
segment_DNAcopy <- function(Log2) {
    ## If segmentation is required, DNAcopy is a good choice. Must be installed. 
    #library(DNAcopy)
    # cat('..Using DNAcopy to create segments:\n')
    segs=NULL
    chroms=c(as.character(1:22),'X','Y')
    chroms=paste('chr',chroms,sep='')
    for (c in 1:24) { # segment chromosome
        tlog=Log2[Log2$Chromosome==chroms[c],]    ## Log-R of this chromosome
        if (nrow(tlog)>0) {                        ## (ChrY may be absent)
            cnaObject=segment(smooth.CNA(CNA(tlog$Value, rep(c,nrow(tlog)), tlog$Start, data.type='logratio',sampleid=paste('chr',c))), undo.splits='sdundo', undo.SD=1,verbose=0) 
            segs=rbind(segs,cnaObject$output)        ## Add result to data frame
        }
    }
    colnames(segs)=c('x','Chromosome','Start','End','Markers','Value')
    segs$Chromosome=chrom_ucsc(segs$Chromosome)    ## Adjust chromosome names to standard
    return(segs[,2:6]) 
}
###
readLog2 <- function() {
    ## This function reads Log-ratio from the file "probes.txt" which must be present in the current directory.
    # Log2=NULL
    # try( Log2 <- read.csv(file='probes.txt',header=T,sep='\t'), silent=T)
    # if (is.null(Log2)) {
    #     try( Log2 <- read.csv(file='_probes.txt',header=T,sep='\t'), silent=T)
    #     if (!is.null(Log2)) cat(' ..found _probes.txt')
    # }

    if(file.exists('probes.txt')) {
        Log2 <- read.csv(file='probes.txt',header=T,sep='\t')
    } else if(file.exists('_probes.txt')) {
        Log2 <- read.csv(file='_probes.txt',header=T,sep='\t')
    } else if(file.exists('raw.txt')) {
        Log2 <- read.csv(file='raw.txt',header=T,sep='\t')
        colnames(Log2)[4:5] <- c('Start','Value')
        Log2$Chromosome <- paste('chr',Log2$Chromosome,sep='')
        Log2$Chromosome[Log2$Chromosome == 'chr24'] <- 'chrX'
        Log2$Chromosome[Log2$Chromosome == 'chr25'] <- 'chrY'
    } else {
        print('No probes.txt found!')
    }

    Log2$Chromosome <- as.character(Log2$Chromosome)
    #c("Chromosome", "Start", "End", "Value", "Array")


    ## This code was used if Log-R must be read from .CNCHP file (Affymetrix Genotyping Console or APT). 
    ## NOT currently supported downstream as .CNCHP is lacks allele-specific information for Affy 250k/500k
    # if (is.null(Log2)) {   
    # dir=dir()
    # ix=grep('cnchp',dir,ignore.case=TRUE)
    # if (length(ix)==1)
    # cat(' ..found',dir[ix])
    # library(affxparser)
    # temp=readCcg(dir[ix])
    # temp=temp$dataGroups$MultiData$dataSets$CopyNumber$table
    # temp$Chromosome[temp$Chromosome==24]=23 # strange numbers in cnchp
    # temp$Chromosome[temp$Chromosome==25]=24
    # Log2 = data.frame('Chromosome'=chrom_ucsc(temp$Chromosome), 'Start'=temp$Position, 'End'=temp$Position, 'Value'=temp$Log2Ratio)
    # save.txt(Log2, file='_probes.txt')
    # cat(' ..wrote _probes.txt')
    # adif=NULL
    # try(
    # adif <- temp$Allele,
    # silent=T
    # )      
    # if (!is.null(adif)) {
    # ix= !is.nan(adif)
    # alf=data.frame('Chromosome'=chrom_ucsc(temp$Chromosome[ix]), 'Start'=temp$Position[ix], 'End'=temp$Position[ix], 'Value'=0.5+adif[ix])
    # save.txt(alf, file='_snps.txt')
    # cat(' ..wrote _snps.txt')
    # }
    # }
    return (Log2) 
}
###
readAlf <- function(localDir=NULL) {
    ## This funciton reads allele frequency [B/(A+B)] from the file 'snps.txt', which must be present in the current directory.
    # alf=NULL
    # try( alf <- read.csv(file='snps.txt',header=T,sep='\t'), silent=T)
    # if (is.null(alf)) {
    #     try( alf <- read.csv(file='_snps.txt',header=T,sep='\t'), silent=T)
    #     if (!is.null(alf)) cat(' ..found _snps.txt')
    # }

    if(file.exists('snps.txt')) {
        alf <- read.csv(file='snps.txt',header=T,sep='\t')
    } else if(file.exists('_snps.txt')) {
        alf <- read.csv(file='_snps.txt',header=T,sep='\t')
    } else if(file.exists('raw_snps.txt')) {
        alf <- read.csv(file='raw_snps.txt',header=T,sep='\t')
        colnames(alf)[c(4,6)] <- c('Start','Value')
        alf$Chromosome <- paste('chr',alf$Chromosome,sep='')
        alf$Chromosome[alf$Chromosome == 'chr24'] <- 'chrX'
        alf$Chromosome[alf$Chromosome == 'chr25'] <- 'chrY'
    } else {
        print('No snps.txt found!')
    }
    #c("Chromosome", "Start", "End", "Value", "Array")

    return (alf) 
}
###
readSegments <- function() {
    ## This function reads segment information in 'segments.txt' which must be present in the current folder.
    ## The author recommends SNP-rank segmentation (NEXUS) or another CBS such as that in DNACopy.
    ## Using a HMM is not recommended unless you have a homogenous, diploid sample. (And then there is more user-friendly software anyway.)
    segments=NULL
    if(file.exists('segments.txt')) {
        segments <- read.csv(file='segments.txt',header=T,sep='\t')
    }
    # try( segments <- read.csv(file='segments.txt',header=T,sep='\t'),silent=T)
    # if (!is.null(segments)) cat(' ..found segments.txt')
    if (is.null(segments)) { 
        # try( segments <- read.csv(file='_segments.txt',header=T,sep='\t'),silent=T)
        # if (!is.null(segments)) cat(' ..found _segments.txt')
    }
    if(!is.null(segments)) {
        segments$Chromosome <- as.character(segments$Chromosome)
    }
    return (segments)
}
###
chrom_n <- function(data) {
    out=rep(Inf,length(data))
    for (c in c(1:22)) {
        out[data==paste('chr',c, sep='')]=c
    }
    out[data=='chrX']=23
    out[data=='chrY']=24
    out[data=='chrM']=25
    return(out)
}
makeRegions <- function(Log2, alf, segments,dataType='Nexus',min=30,matched=FALSE,allelePeaks=F) {
    ## makeRegions is similar to "regsfromsegs" except regions are not subdivided before calculation of mean Log-R and Allelic Imbalance Ratio.
    #browser()
    regions=segments
    regions$Chromosome=as.character(segments$Chromosome)            ## Chromosome
    regions$lengthMB=round((regions$End-regions$Start)/1000000,3)    ## length in megabases
    regions$probes=0                                                ## # of probes
    regions$snps=0                                                    ## # of bi-allelic probes (SNPs)
    #regions$log2=round(regions$Value,4)
    regions$imba=NA                                                 ## Allelic Imbalance Ratio
    regionIx=NULL                                                   ## Not currently used
    regionIx$Log2 <- list()
    regionIx$alf <- list()
    
    pMatch=segMatch(regions,Log2)
    na=is.na(pMatch[,1])
    regions=regions[!na,]
    pMatch=pMatch[!na,]
    sMatch=segMatch(regions,alf)
    
    for (i in 1:nrows(regions)) {
        #log2temp=which(equals(Log2$Chromosome,regions$Chromosome[i])) ## index of Log-R (current chrom)
        #alftemp=which(equals(alf$Chromosome,regions$Chromosome[i]))    ## index of Allele frequency (current chrom)
        
        #log2temp=log2temp [Log2$Start[log2temp]>=regions$Start[i] & Log2$Start[log2temp]<regions$End[i]]    ## index of Log-R (current segment)
        #alftemp=alftemp [alf$Start[alftemp]>=regions$Start[i] & alf$Start[alftemp]<regions$End[i]]        ## index of Allele frequency (current segment)
        
        log2temp=Log2$Value[pMatch[i,1]:pMatch[i,2]]    ## index of Log-R (current segment)
        alftemp=numeric(0)
        if (!is.na(sMatch[i,1])) alftemp=alf$Value[sMatch[i,1]:sMatch[i,2]]        ## index of Allele frequency (current segment)
        
        regions$probes[i]=length(log2temp)
        regions$snps[i]=length(alftemp)
        #regionIx$Log2[[i]]=log2temp            ## indexes of Log-R and Allele Frequency are saved for future use (plot color coding)
        #regionIx$alf[[i]]=alftemp
        #log2temp=Log2$Value[log2temp]
        regions$log2[i]=median(log2temp, na.rm=T)
        #alftemp=alf$Value[alftemp]
        temp = NA
        try(temp <- allelicImbalance(alftemp,min,matched=matched,allelePeaks=allelePeaks),silent=T)
        regions$imba[i]  <- temp

        # if (length(alftemp)>3) {                ## Prepare to calculate Allelic Imbalance Ratio
        #     if (dataType=="Nexus") t1=sort( abs(alftemp-0.5) ) ## distance from middle (het), ascending
        #     if (dataType=="CNCHP") t1=sort( abs(alftemp)-0 ) ## This is not currently in use, was intended for analysis of SNP6 .CNCHP data.
        #     if (length(unique(t1))>5) { ## Avoid calculating Allelic Imbalance Ratio unless there are several different values to cluster
        #         xx=NA
        #         xx=try(kmeans(t1, 2), silent=T)            ## This part is nearly identical to that of 'regsfromsegs()'
        #         try( if (min(xx$size) > 0.05*max(xx$size)) {
        #             xx=xx$centers
        #         } else xx=NA, silent=T)    
        #     } else xx=NA
        #     try(regions$imba[i] <- round( min(xx)/max(xx) ,2), silent=T)
        # }
        # if(matched == T | is.na(regions$imba[i])) {
        #     regions$imba[i] <- 2*median(abs(alftemp-.5),na.rm=T)
        # }
    }
    return(list('regions'=regions,'regionIx'=regionIx))
}
###
equals <- function(a,b) {
    ## a helper function in case factors are compared
    a <- as.character(a)
    b <- as.character(b)
    return (a==b)
}
###
nrows <- function(data) dim(data)[1] ## don't ask me about this one.
###
load.txt <- function(file) { ## just to be rid of "sep"
    read.csv(file,sep='\t') 
}
###
save.txt <- function(data,file='temp', row.names=F, col.names=T, append=F) {   ## simplified for convenience
    write.table(data,file=file,sep='\t',quote=F,row.names=row.names, col.names=col.names, append=append)
}
###
deChrom_ucsc <- function(data) {
    ##  chroms in 1:24 and chr1:chrY format conversion
    if (length(data)==0) return (data)
    keep_index=NULL
    CHR=rep(0,length(data))
    existing=unique(data)
    for (i in 1:length(existing))   {
        temp=strsplit(as.character(existing[i]),'_')
        temp=strsplit(temp[[1]],'chr')[[1]][2]
        if (temp=='X') temp='23'
        if (temp=='_X') temp='23'
        if (temp=='Y') temp='24'
        CHR[data==existing[i]]=as.integer(temp)
    }
    return (CHR)
}
###
chrom_ucsc <- function(data) {
    ## same as above, backwards
    n=length(data)
    out=rep("",n)
    for (i in 1:n) {
        temp=as.character(data[i])
        if (temp=='23') temp='X'
        if (temp=='24') temp='Y'
        out[i]=paste('chr',temp, sep='')
    }
    return(out)
}
###
subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y) ## for subplots
vplayout <- function(x, y) {        
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(y,x)))
}
###
getSubdirs <- function() {
    ## Reach for the subdirectories of current directory.
    dir=dir()
    here=getwd()
    n=length(dir)
    subdir=F
    dir_ix=rep(FALSE,n)
    for (i in 1:n) {
        try(setwd(dir[i]), silent=T)
        if (getwd()!=here) { # in another directory
            setwd('..')
            subdir=T
            dir_ix[i]=T
        }
    }
    if (subdir) return (dir[dir_ix])
    return (NULL)
}
###
thisSubdir <- function() { 
    ## get the name of this subdirectory
    here=getwd()
    here=strsplit(here,'/')[[1]]
    here=here[length(here)]
    return(here)
}
###
mySorter <- function(data) {
    ## sort data on $Chromosome and $Start
    temp=data[order(data$Start),]
    t=deChrom_ucsc(temp$Chromosome)
    temp=temp[!is.na(t),]
    t=t[!is.na(t)]
    temp=temp[order(t),]
    return (temp)
}
###
probesInRegions <- function(regions,probes) {
    ## This function is not in use.
    rchr=as.character(regions$Chromosome)
    pchr=as.character(probes$Chr)
    nprobes=0
    
    for (i in 1:length(rchr)) {
        t=strsplit(rchr[i],'chr')[[1]][2]
        here=probes[(t==pchr),] #finding probe subset on this chr
        matching=(here$Position>=regions$Start[i]) & (here$Position<=regions$End[i])
        nprobes[i]=sum(matching)
    }
    return(nprobes)
}
###
getOverlap <- function(s1,e1,s2,e2) { ## This is not currently used anywhere - it is destined for group studies and pileups.
    return(max((min(e1,e2)-max(s1,s2)),0))
}
###
weightedMedian <- function(data,weights) {
    try( if (length(weights)==0) return (NULL), silent=T)
    try( return (median(rep(data,weights),na.rm=T)), silent=T)
    return (NULL)
}
weightedMean <- function(data,weights) {
    if (length(weights)==0) return (NULL)
    return (mean(rep(data,weights),na.rm=T))
}
###
myMeans <- function(data,k) { ## Not currently used
    best <- 0
    clus <- list()
    for (i in 1:10) {
        clus[[1]] <- kmeans(data,k)
        var <- var(clus[[1]]$centers)
        if (var>best) {
            ix <- i
            best <- var
        }
    }
    return (clus[[ix]])
}
is.even <- function(data) {
    return (trunc(data/2)*2 == data)
}
is.odd <- function(data) {
    return (trunc(data/2)*2 != data)
}
neighbour <- function(data) {
    n <- length(data)
    dif <- data[2:n]-data[1:(n-1)]
    return (mean(dif))
}
is.autosome <- function(vector) {
    temp <- deChrom_ucsc(vector)
    return (temp<23)
}

## 
###
findCNs <- function(Log2,alf,allRegions,regs,name=thisSubdir(),maxCn=10,ceiling=1,sampleInfo=NULL) {
    ## This function takes an estimate of the Log-R of copy numbers 1, 2 and 3. At least two of these should be entered.
    ## Then, the Log-R and Allelic Imbalance Ratio of all possible copy number variants up to maxCn are estimated from
    ## the Log-R and Allelic Imbalance Ratio of all the segments. 
    if (is.null(sampleInfo)) cat ('there was no estimation available for',name)
    
    cns=1:maxCn; est=sampleInfo[1:3]; est[est==' ']=NA; est=as.numeric(as.character(est))
    m <- lm(2^est ~ cns[1:3])$coefficients # can handle one NA in est. This model is for "ratio as a function of copy number".
    est[is.na(est)] = log2(m[1]+cns[which(is.na(est))]*m[2]) # simple linear regression to fill the missing
    
    tix=NULL     #temporary index
    probes=NULL  #number of probes at each copy number, for weighting
    int=NULL     ## contains Log-R estimate of each (total) copy number
    ai=NULL         ## contains Allelic Imbalance Ratio estimate of each copy number variant.
    regions <- allRegions$regions
    regions <- regions[(is.autosome(regions$Chromosome)&regions$lengthMB>1)&(!is.na(regions$imba)),] ## will use these regions
    
    ## likely cn2 regions sit near the estimate.
    expectedAt <- est[2]
    tix$cn2 <- abs(regions$log2 - expectedAt) < diff(est)[2]/3    ## index of likely cn2 regions
    temp <- regions[tix$cn2,]                                ## cn2 regions
    med <- weightedMedian(temp$log2,temp$probes)            ## improved value of Log-R at cn2 (returns NULL if theres nothing there)
    probes[2] <- sum(temp$probes)
    int[2] <- ifelse(!is.null(med),med,expectedAt)            ## saved to int.
    
    ## likely cn1 regions sit near estimate                                         
    expectedAt <- est[1]                                     ## cn1 is expected here
    tix$cn1 <- abs(regions$log2 - expectedAt) < diff(est)[1]/3        ## index of likely cn1 regions
    temp <- regions[tix$cn1,]
    med <- weightedMedian(temp$log2,temp$probes)
    probes[1] <- ifelse(is.null(med), 50, sum(temp$probes))
    int[1] <- ifelse(!is.null(med),med,expectedAt)
    
    ## likely cn3 regions sit near estimate
    expectedAt <- est[3]
    tix$cn3 <- abs(regions$log2 - expectedAt) < diff(est)[2]/3
    temp <- regions[tix$cn3,]
    med <- weightedMedian(temp$log2,temp$probes)
    probes[3] <- ifelse(is.null(med), 50, sum(temp$probes))
    int[3] <- ifelse(!is.null(med),med,expectedAt)
    
    ## cn4 follows at ...
    m <- lm(2^int[1:3] ~ cns[1:3], weights=probes[1:3])$coefficients ## use regression to estimate.
    expectedAt <- log2(m[1]+4*m[2])
    tix$cn4 <- abs(regions$log2 - expectedAt) < mean(diff(int))/3
    temp <- regions[tix$cn4,]
    med <- weightedMedian(temp$log2,temp$probes)
    probes[4] <- ifelse(is.null(med), 50, sum(temp$probes))
    int[4] <- ifelse(!is.null(med),med,expectedAt)
    
    ## generalized for higher cns
    for (cn in 5:maxCn) {
        thisCn <- paste('cn',cn,sep='')
        prevCn <- paste('cn',cn-1,sep='')
        pprevCn <- paste('cn',cn-2,sep='')
        m <- lm(2^int[1:(cn-1)] ~ cns[1:(cn-1)], weights=probes[1:(cn-1)])$coefficients
        expectedAt <- log2(m[1]+cn*m[2])
        tix[[thisCn]] <- abs(regions$log2 - expectedAt) < mean(diff(int))/5
        temp <- regions[tix[thisCn][[1]],]
        med <- weightedMedian(temp$log2,temp$probes)
        probes[cn] <- ifelse(is.null(med), 0, sum(temp$probes))
        int[cn] <- ifelse(!is.null(med),med,expectedAt)
    }
    
    ## likely cn0 regions sit below cn1 - delta:
    expectedAt <- log2(m[1]+0*m[2])
    tix$cn0 <- abs(regions$log2 - expectedAt) < 0.5*(int[2]-int[1])
    temp <- regions[tix$cn0,]
    med <- weightedMedian(temp$log2,temp$probes)
    int0 <- ifelse(!is.null(med),med,expectedAt) ## "int0"
    
    
    ## Estimate tumor dna content from intensity-cn relationship and average ploidy --UNRELIABLE
    md <- lm(2^int ~ cns, weights=probes)$coefficients
    m=NULL; m$intercept=md[1][[1]]; m$k=md[2][[1]]
    probes[is.na(probes)]=0
    m$meanCn <- mean(rep(cns, probes), na.rm=T)
    m$theoretical_delta=1/m$meanCn
    #m$real_delta=0.57*m$theoretical_delta ## The 0.57 is empirical from cancer cell lines    
    #m$dnafrac=m$k/m$real_delta
    #m$cellfrac=1/(1+m$meanCn/2*(1/m$dnafrac-1))    
    

    
    loh_exp <- as.numeric(sampleInfo[4])
    ## at cn2, find the variant clusters (normal and CNNLOH)
    ix <- (abs(regions$log2 - int[2]) < 0.2*(int[3]-int[2]) ) # taking only closely-matching segments
    data <- regions[ix,]
    data <- data[!is.na(data$imba),] # ...with a calculated allelic imbalance.
    ## 2m1
    expectedAt=0.1
    ix <- abs(data$imba-.1) < abs(data$imba-loh_exp)
    med <- weightedMedian(data$imba[ix],data$snps[ix]) ## Average allelic imbalance weighted on snp count
    ai$cn2m1 <- ifelse (!is.null(med),med,expectedAt)
    ## loh
    expectedAt=loh_exp
    ix <- !ix
    med <- weightedMedian(data$imba[ix],data$snps[ix]) 
    ai$cn2m0 <- ifelse (!is.null(med),med,expectedAt)
    
    
    ## for cn1 (and 0)
    ix <- (abs(regions$log2 - int[1]) < 0.2*(int[2]-int[1]) )
    data <- regions[ix,]
    data <- data[!is.na(data$imba),]
    expectedAt <- (ai$cn2m0+ai$cn2m1)*3/5     ## Decent estimate.
    med <- weightedMedian(data$imba,data$snps) ## Average allelic imbalance weighted on snp count
    ai$cn1m0 <- ifelse (!is.null(med),med,expectedAt)    ## Will be NA if there was no CNNLOH
    ai$cn0m0 <- NA #unimportant
    
    ## for cn3:
    ix <- (abs(regions$log2 - int[3]) < 0.2*(int[4]-int[3]) )
    data <- regions[ix,]
    data <- data[!is.na(data$imba),]
    
    if (!is.na(ai$cn2m0)) { # CNNLOH is not missing
        range <- ai$cn2m0-ai$cn2m1 #the distance between 2normal and CNNLOH
        # get the 3(1) regions: 
        expectedAt <- ai$cn2m1+range/3 # this is an approx if cn3m1 is absent.
        
        ix <- (data$imba<ai$cn2m0) & (data$imba>ai$cn2m1) # take regions with less AI than cn2m0 but more than cn2m1
        med <- weightedMedian(data$imba[ix],data$snps[ix]) # average allelic imbalance weighted on snp count
        ai$cn3m1 <- ifelse (!is.null(med),med,expectedAt)
        
        # now for cn3m0
        expectedAt <- ai$cn3m1 + range*0.9 # approx of cn3m0
        
        ix <- (abs(data$imba-expectedAt) / abs(data$imba-ai$cn3m1)) < 0.5  # take those much closer to exp than to cn3m1
        med <- weightedMedian(data$imba[ix],data$snps[ix]) 
        try (if (med<ai$cn2m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
        ai$cn3m0 <- ifelse (!is.null(med),med,expectedAt)
    } else { ## (if CNNLOH was missing, use a different strategy.)
        data <- data[data$imba>ai$cn2m1,] # avoid 2(0) / 4(0) heterogeneity
        data <- data[data$imba<0.85,] # avoid highly-hom (clus-errors)
        #centers <- sort(kmeans(rep(data$imba,data$snps),2)$centers)## bad idea if either is missing! (common on low-mut samples)
        ix <- data[data$imba<0.4,]
        med <- weightedMedian(data$imba[ix],data$snps[ix])
        ai$cn3m1 <- ifelse (!is.null(med),med,ai$cn2m1+0.1)
        ix <- data[data$imba>0.4,]
        med <- weightedMedian(data$imba[ix],data$snps[ix])
        ai$cn3m0 <- ifelse (!is.null(med),med,ai$cn2m1+0.35)
        ai$cn2m0 <- ai$cn2m1+(ai$cn3m0-ai$cn3m1) # the CNNLOH estimate from CN3
        #print(cat('Warning: cn2m0 missing, estimated cn3m1,cn3m0 centers at',ai$cn3m1,ai$cn3m0,'\n'))
    }
    
    if (is.na(ai$cn1m0)) ai$cn1m0 <- (ai$cn3m1+ai$cn2m0)/2 ## If deletions were missing, place an estimate from cn3 
    
    ## now for cn4
    ix <- abs(regions$log2 - int[4]) < 0.2*(int[4]-int[3]) 
    data <- regions[ix,]
    data <- data[!is.na(data$imba),]
    
    # get the 4(2) regions: they are at AI about cn2m1
    expectedAt <- ai$cn2m1 # this is a good approx 
    
    ix <- data$imba<ai$cn3m1 # let all below cn3m1 in
    med <- weightedMedian(data$imba[ix],data$snps[ix]) 
    ai$cn4m2 <- ifelse (!is.null(med),med,expectedAt)
    if (ai$cn4m2 > ai$cn2m1) ai$cn4m2 <- ai$cn2m1 # don't let it sneak off.'
    
    data <- data[!ix,] # coutinue with the remaining
    # now 4(1) has less ai than 3(0) -> sits at about cn3m1+(cn3m1-cn2m1).
    expectedAt <- ai$cn3m1+(ai$cn3m1-ai$cn2m1)
    ix <- abs(data$imba-expectedAt)<abs(data$imba-ai$cn3m0) # take those closer to exp than to 3(0)
    med <- weightedMedian(data$imba[ix],data$snps[ix]) 
    ai$cn4m1 <- ifelse (!is.null(med),med,expectedAt)
    
    # now 4(0) is at about 3(0) + [0.9 - 3(0)] * [3(0)-2(0)] / [0.9-2(0)] 
    expectedAt <- ai$cn3m0 + (ceiling-ai$cn3m0)*(ai$cn3m0-ai$cn2m0)/(ceiling-ai$cn2m0)
    ix <- abs(data$imba-expectedAt) < 0.2*(expectedAt-ai$cn4m2) # we take those very close to exp
    med <- weightedMedian(data$imba[ix],data$snps[ix])
    try (if (med<ai$cn3m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
    ai$cn4m0 <- ifelse (!is.null(med),med,expectedAt)
    
    ## generalization for higher copy numbers. it gets complicated.
    for (cn in 5:maxCn) {
        thisCn <- paste('cn',cn,sep='')
        prevCn <- paste('cn',cn-1,sep='')
        pprevCn <- paste('cn',cn-2,sep='')
        
        ix <- (abs(regions$log2 - int[cn]) < 0.2*(int[cn]-int[cn-1]) )
        data <- regions[ix,]
        data <- data[!is.na(data$imba),]
        data <- data[data$lengthMB>3,] # long regions for safety
        
        ## try to find variants, starting with LOH
        # LOH such as 5(0)
        mi <- 0
        thisVariant=paste(thisCn,'m',0,sep='')
        c4m0 <- ai[paste(prevCn,'m',mi,sep='')][[1]] # relative naming for clarity
        c3m0 <- ai[paste(pprevCn,'m',mi,sep='')][[1]]
        
        expectedAt <- ceiling-((ceiling-c4m0)*(ceiling-max(c4m0,c3m0))/(ceiling-min(c3m0,c4m0)))
        #ai[thisVariant] <- expectedAt
        # we then take those close to exp (within delta[3(0),4(0)])
        ix <- abs(data$imba-expectedAt) < abs(c4m0 - c3m0) 
        med <- weightedMedian(data$imba[ix],data$snps[ix]) 
        try (if (med<ai$cn4m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
        
        ## then from balanced to less balanced
        minorVariants=trunc(cn/2):1 
        first <- T
        for (mi in minorVariants) {
            thisVariant=paste(thisCn,'m',mi,sep='')
            if (mi==cn/2) {
                # We have balanced variant, so rather easy.
                expectedAt <- ai[paste(pprevCn,'m',mi-1,sep='')][[1]] # this is a good approx, balanced at cn-2
                ix <- data$imba<ai[paste(prevCn,'m',mi-1,sep='')][[1]] # let all below (cn-1, mcn-1) in
                med <- weightedMedian(data$imba[ix],data$snps[ix]) 
                ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
                if (ai[thisVariant] > ai[paste(pprevCn,'m',mi-1,sep='')][[1]]) ai[thisVariant] <- ai[paste(pprevCn,'m',mi-1,sep='')][[1]] # don't let it sneak off.'
                first <- F
            } else if (first) { 
                # its not balanced but its the most balanced of the unbalanced. something like 5(2)
                expectedAt <- 0.5*( ai[paste(prevCn,'m',mi,sep='')][[1]] + ai[paste(pprevCn,'m',mi-1,sep='')][[1]] ) # that means between 4(2) and 3(1)
                ix <- abs(data$imba-expectedAt) < ( ai[paste(prevCn,'m',mi,sep='')][[1]] - ai[paste(pprevCn,'m',mi-1,sep='')][[1]] ) /3 # let all "between" 4(2) and 3(1) in
                med <- weightedMedian(data$imba[ix],data$snps[ix]) 
                ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
                first <- F
            } else {
                # not the most balanced unbalanced variant, for example 5(1):
                expectedAt <- ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]] + 
                    (minorVariants[1]-mi) * (ai[paste(thisCn,'m',0,sep='')][[1]] - 
                                                 ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]]) / trunc(cn/2) # 5(2) + (which)* 5(0)-5(2) /(n)
                ix <- abs(data$imba-expectedAt) < ( ai[paste(prevCn,'m',mi,sep='')][[1]] - ai[paste(pprevCn,'m',mi-1,sep='')][[1]] ) /3 # let all "between" 4(1) and 3(0) in
                med <- weightedMedian(data$imba[ix],data$snps[ix]) 
                ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
            } 
        } # done with minor variants
    } # done with copy numbers
    int=as.list(c(int0,int)); names(int)=paste('cn',0:maxCn,sep='')
    return(list('int'=int,'ai'=ai,'model'=m))
}



###
setCNs <- function(allRegions,int,ai,model,maxCn=12) {
    ##  Assign total and minor copy numbers to all segments.
    regions <- allRegions$regions[,-4]    ## This time, work on all segments available.
    
    Cn <- NULL            ## Total copy number
    mCn <- NULL            ## Minor allele copy number
    Cnx <- NULL        ## Variant label. ('cnXmY')
    
    intDist <- NULL        ## distance to certain Log-R
    imbaDist <- NULL    ## distance to certain allelic imbalance
    
    ## give each segment the best matching CN and mCN
    for (i in 1:nrow(regions)) {
        # set total copy number
        distance <- Inf
        for (cn in 0:maxCn) {
            t_int <- int[[paste('cn',cn,sep='')]]   ## get Log-R of particular cn from 'int'
            t_dis <- abs(regions$log2[i]-t_int)            ## distance to that particular cn
            try (if (t_dis < distance) {                        ## nearest so far, save.
                distance <- t_dis -> intDist[i]
                Cn[i] <- cn
            }, silent=T)
        }
    }        

    ### Calculate model based Cns
    Cnx=Cn; for (cn in 0:maxCn) {
        ix <- Cn==cn
        ix[is.na(ix)] <- F
        Cnx[ix]=(cn + (2^regions$log2[ix]-2^int[[paste('cn',cn,sep='')]])/model$k)
    }; Cnx[Cnx<0]=0
    Cn=round(Cnx)
    
    ## Set minor CN
    for (i in (1:nrow(regions))[!is.na(Cn)]) {
        # set minor CN
        distance <- Inf
        if (Cn[i]<=1) {
            mCn[i] <- 0
        } else if (Cn[i]<=maxCn & !is.na(regions$imba[i])) for (m in 0:trunc(Cn[i]/2)) {
            t_ai <- ai[paste('cn',Cn[i],'m',m,sep='')][[1]]
            t_dis <- abs(regions$imba[i]-t_ai)
            #cat('Line',i,'Cn', Cn[i],'m:',m,'Comparing',t_dis,'with',t_int,'.......')
            try( if (t_dis < distance) {
                distance <- t_dis -> imbaDist[i]
                mCn[i] <- m
            }, silent=T)
        } else mCn[i] <- NA
        #fullCN[i] <- paste('cn',Cn[i],'m',mCn[i],sep='') # Full description
    }
    ### Do not calculate mCn for high copy numbers
    mCn[Cn>maxCn] <- NA
    
    
    ## Special treatment for the X chromosome if male sample: work in progress
    #temp=regions[regions$Chromosome=='chrX' & regions$lengthMB>5,]
    
    
    regions$Cn <- Cn
    regions$mCn <- mCn
    regions$Cn_ <- round(Cnx,2)
    #regions$fullCN <- fullCN

    ## Merge consecutive identical variants
    data <- regions[,c(-7,-8,-11)]
    fullCN=paste(data$Cn,data$mCn,sep='m')
    row <- 1
    while (row < nrow(data)) { # while not on the last row
        #print(row)
        if (((data$Chromosome[row] == data$Chromosome[row+1]) & (fullCN[row] == fullCN[row+1])) & ( (data$Start[row+1]-data$End[row])<5000 ) ) { 
            ## segments are adjacent with same copy number and not separated by 5kb+ (centromere)
            data$End[row] <- data$End[row+1]
            data$probes[row] <- data$probes[row]+data$probes[row+1]
            data$snps[row] <- data$snps[row]+data$snps[row+1]
            data$lengthMB[row] <- data$lengthMB[row]+data$lengthMB[row+1]
            data <- data[-(row+1),]  # Merges and deletes row
        } else {
            row <- row+1 # or leaves it
        }
    }
    
    ### beräkna tumorcellshalt!!
    # Mha modellparameter
    ## empirical delta-R for ±one copy in diploid sample: 0.63 for Cyto/SNP6 and 0.48 for Illu
    model$k2=model$k*model$meanCn
    model$tumorDNA_fromLogR=model$k2/0.63
    model$tumorCell_fromLogR=1/ (1+model$meanCn/2* (1/model$tumorDNA_fromLogR -1))
    
    # Mha SNPs
    cn_near=c(floor(model$meanCn),ceiling(model$meanCn))
    ai_near=c(ai[[paste('cn',cn_near[1],'m0',sep='')]],ai[[paste('cn',cn_near[2],'m0',sep='')]])
    base_ai = ai_near[1] + (model$meanCn-trunc(model$meanCn))*diff(ai_near)
    
    #---The reference data (H1395):
    ref_ai=c(.8,.5,.33,.2,.1)
    ref_TumorCell=c(1,.7,.5,.3,0)
    ref_meanCn=2.5
    ref_tumorDNA=ref_TumorCell*ref_meanCn/(ref_TumorCell*ref_meanCn + (1-ref_TumorCell)*2)
    m2=lm(ref_tumorDNA ~ poly(ref_ai, 2, raw=TRUE))
    #---
    
    model$tumorDNA_fromAlf <- m2[[1]][[1]] + base_ai*m2[[1]][[2]] + base_ai^2*m2[[1]][[3]]
    model$tumorCell_fromAlf=1/ (1+model$meanCn/2* (1/model$tumorDNA_fromAlf -1))
    
    return (list('regions'=regions,'merged'=data,'model'=model))
}


karyotype_check <- function(chr,start,end,int,ai,Cn,mCn,t,name='') { #xlim=c(-1.02,1.02),ylim=0:1) {
    ## TAPS scatter plot of a full sample, used for visual quality control. 
    
    png(paste(name,'.karyotype_check.png',sep=''),width=1300,height=1300)
    
    #ideogram=getIdeogram()
    colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")
    
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
    
    
    
    variants <- unique(labels)                ## The variants present in this data
    variants_data <- NULL    
    for (v in variants) {                            ## Extracts the Log-R and Allelic Imbalance Ratio of variants
        t_cn <- paste('cn',strsplit(v,'m')[[1]][1],sep='')
        t_int <- as.numeric(t$int[t_cn][[1]])
        t_ai <- as.numeric(t$ai[paste('cn',v,sep='')][[1]])
        variants_data <- rbind(variants_data,c(t_int,t_ai))
    }
    
    
    col <- rep('#B0B0B070',length(chr))
    
    plot(int,ai,
         pch=16,
         cex=c(size),
         cex.lab=2,
         mar=c(0.1,0.1,0.1,0.1),
         main = "",
         xlab = "",
         ylab = "",
         col = col,
         xlim=c(-1,2),ylim=c(0,1))
    
    text(variants_data[,1],variants_data[,2],
         labels=variants,
         col='black',
         cex=2
    )
    
    mtext("Allelic imbalance",side=2,line=2,cex=2)
    mtext("Average log-ratio",side=1,line=2.3,cex=2)
    mtext(name,side=3,line=1.5,cex=2)


    dev.off()
    
}

# getIdeogram <- function() {
#     c=1:24
#     chr=as.character(c)
#     chr[23:24]=c('X','Y')
#     chr=paste('chr',chr,sep='')
#     length=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 
#              146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 
#              102531392,  90354753,  81195210,  78077248,  59128983,  63025520,  48129895,  
#              51304566, 155270560,  59373566)
#     start=c(121500000,  90500000,  87900000,  48200000,  46100000,  58700000,  58000000,  
#             43100000,  47300000,  38000000,  51600000,  33300000,  16300000,  16100000,  
#             15800000,  34600000,  22200000,  15400000,  24400000,  25600000,  10900000,  
#             12200000,  58100000,  11600000)
#     mid=c(1.25e+08, 9.33e+07, 9.10e+07, 5.04e+07, 4.84e+07, 6.10e+07, 5.99e+07, 4.56e+07, 
#           4.90e+07, 4.02e+07, 5.37e+07, 3.58e+07, 1.79e+07, 1.76e+07, 1.90e+07, 3.66e+07, 
#           2.40e+07, 1.72e+07, 2.65e+07, 2.75e+07, 1.32e+07, 1.47e+07, 6.06e+07, 1.25e+07)
#     end=c(128900000,  96800000,  93900000,  52700000,  50700000,  63300000,  61700000,  
#           48100000,  50700000,  42300000,  55700000,  38200000,  19500000,  19100000,  
#           20700000,  38600000,  25800000,  19000000,  28600000,  29400000,  14300000,  
#           17900000,  63000000,  13400000)
#     return(data.frame(c,chr,length,start,mid,end))
# }





## Pileups
pileup <- function(regions) { #, nsamples) {
    if (nrow(regions)==0) return(data.frame(Chromosome=integer(0),Start=integer(0),End=integer(0),Count=integer(0)))
    chrs <- unique(as.character(regions$Chromosome))
    if (length(chrs)>1) { # if more than one chrom, run them separately.
        pile=NULL
        for (i in 1:length(chrs)) pile <- rbind(pile,pileup(regions[as.character(regions$Chromosome)==chrs[i],]))
        #pile$Percent <- round(100*pile$Count/nsamples) # this is done once.
        return(pile[order(pile$Chromosome),])
    } else {
        c=as.character(regions$Chromosome[1])
        starts=data.frame(Chromosome=c,Start=regions$Start,End=regions$Start,type='start')
        ends=data.frame(Chromosome=c,Start=regions$End,End=regions$End,type='end')
        pile <- rbind(starts,ends)
        pile <- pile[order(pile$Start),]
        pile$Count <- 1
        
        for (i in 2:nrow(pile)) {
            pile$Count[i] <- pile$Count[i-1] + ifelse(pile$type[i]=='start',1,-1)
            pile$End[i-1] <- pile$Start[i]
        }
        pile <- pile[(pile$End > pile$Start) & (pile$Count > 0),]
        return(pile)
    }
}

pileup_dif <- function(regions1,regions2) {
    chrs <- unique(c(as.character(regions1$Chr),as.character(regions2$Chr)))
    if (length(chrs)==0) return(data.frame(Chromosome=integer(0),Start=integer(0),End=integer(0),Count1=integer(0),Count2=integer(0)))
    if (length(chrs)>1) { # if more than one chrom, run them separately.
        pile=NULL
        for (i in 1:length(chrs)) pile <- rbind(pile,pileup_dif(regions1[as.character(regions1$Chr)==chrs[i],], regions2[as.character(regions2$Chr)==chrs[i],]))
        #pile$Percent <- round(100*pile$Count1/n1 - 100*pile$Count2/n2) # this is done once.
        #pile <- pile[pile$Percent != 0,]
        return(pile[order(pile$Chr),])
    } else {
        c <- chrs
        starts1 <- starts2 <- ends1 <- ends2 <- NULL
        if (nrow(regions1)>0) {
            starts1=data.frame(Chromosome=c,Start=regions1$Start,End=regions1$Start,type='start',grp=1)
            ends1=data.frame(Chromosome=c,Start=regions1$End,End=regions1$End,type='end',grp=1)
        }
        if (nrow(regions2)>0) {
            starts2=data.frame(Chromosome=c,Start=regions2$Start,End=regions2$Start,type='start',grp=2)
            ends2=data.frame(Chromosome=c,Start=regions2$End,End=regions2$End,type='end',grp=2)
        }
        pile <- rbind(starts1,ends1,starts2,ends2)
        pile <- pile[order(pile$Start),]
        pile$Count1 <- ifelse(pile$grp[1]==1,1,0)
        pile$Count2 <- ifelse(pile$grp[1]==2,1,0)
        
        for (i in 2:nrow(pile)) {
            if (pile$type[i]=='start') {
                pile$Count1[i] <- pile$Count1[i-1] + ifelse(pile$grp[i]==1,1,0)
                pile$Count2[i] <- pile$Count2[i-1] + ifelse(pile$grp[i]==2,1,0)
            } else {
                pile$Count1[i] <- pile$Count1[i-1] - ifelse(pile$grp[i]==1,1,0)
                pile$Count2[i] <- pile$Count2[i-1] - ifelse(pile$grp[i]==2,1,0)
            }
            pile$End[i-1] <- pile$Start[i]
        }
        
        pile <- pile[(pile$End > pile$Start),]
        return(pile)
    }
}


pvals <- function(regions,tot1,tot2) {
    if (nrow(regions)==0) return(data.frame(Chromosome=integer(0),Start=integer(0),End=integer(0),Count1=integer(0),Count2=integer(0),
                                            percent1=integer(0),precent2=integer(0),precent=integer(0),p.value=integer(0),
                                            odds.ratio=integer(0),conf_low=integer(0),conf_high=integer(0)))
    p <- rep(1,nrow(regions))
    or <- rep(1,nrow(regions))
    conf_low <- rep(NA,nrow(regions))
    conf_high <- rep(NA,nrow(regions))
    for (i in 1:nrow(regions)) {
        temp <- matrix(c(regions$Count1[i],regions$Count2[i],tot1-regions$Count1[i],tot2-regions$Count2[i]),
                       nrow=2,
                       dimnames=list(c('rec','ej'),c('have','nohave')))
        temp <- fisher.test(temp)
        p[i] <- temp$p.value
        or[i] <- temp$estimate
        conf_low[i] <- temp$conf.int[1]
        conf_high[i] <- temp$conf.int[2]
    }
    regions$percent1 <- round(100*regions$Count1/tot1)
    regions$percent2 <- round(100*regions$Count2/tot2)
    regions$percent <- regions$percent1-regions$percent2
    regions$p.value <- p 
    regions$odds.ratio <- or
    regions$conf_low <- conf_low
    regions$conf_high <- conf_high
    return (regions)#[order(regions$p.value,decreasing=F),])
}

getGenes <- function(regs,genes) {
    n <- nrow(regs)
    res <- rep('NA',n)
    for (i in 1:n) {
        ix <- genes$Chromosome==regs$Chromosome[i]
        ix2 <- (genes$End>regs$Start[i])&(genes$Start<regs$End[i])
        ix <- ix & ix2
        temp <- genes$Name[ix]
        if (length(temp)==0) temp <- '-'
        temp <- paste(unique(temp),collapse=',')
        res[i] <- temp
    }
    return(res)
}

vectorMatch <- function(v1,v2) {
    v1 <- as.character(v1); v2 <- as.character(v2)
    n1 <- length(v1); n2 <- length(v2) ## vector lengths
    
    ix1 <- order(v1); ix2 <- order(v2) ## how to sort them
    ori1 <- (1:n1)[ix1]; ori2 <- (1:n2)[ix2] ## how to restore the sorted
    
    v1 <- v1[ix1]; v2 <- v2[ix2] ## now sort
    mix <- rep(NA,n1) ## for each in sorted1, which is it in original 2?
    i1=1; i2=1
    while((i1 <= n1) & (i2 <= n2)) {
        if (v1[i1]>v2[i2]) {
            i2 <- i2+1
        } else if (v1[i1]<v2[i2]) {
            i1 <- i1+1
        } else { ## we have a match. 
            mix[i1] <- ori2[i2]       ## mix[i1] must point to the original pos of i2, in vector 2.
            i1 <- i1+1
            i2 <- i2+1
        }
    }
    return (mix[order(ori1)])
}

pileup_pos <- function(chr,pos,names,chroms, binsize=5000000) {
    data1 <- data.frame(chr=as.character(chr),pos=as.integer(pos),name=names)
    data <- unique(data1[,1:2])
    
    allbins <- NULL
    #binsize <- 5000000 ###warning
    for (i in 1:24) {
        chr <- as.character(chroms$Chr[i])
        pos <- data$pos[data$chr==chr]
        breaks <- c(0,seq(0,chroms$length[i],binsize)+binsize)
        histo <- hist(pos,breaks=breaks,plot=F)
        start <- breaks[-length(histo$breaks)]
        end <- breaks[-1]
        mid=histo$mid
        count <- histo$count
        
        # now add number of samples that make up the count (UNFINISHED)
        for (j in 1:length(start)) { # for every segment
            ix <- pos>=start[j] & pos<=end[j] # which positions are in
            samples[j] <- length(unique(names[ix]))
        }
        
        allbins <- rbind(allbins,data.frame(chr,start,mid,end,count))
    }
    e <- 0 -> s
    for (i in 1:24) s[allbins$chr==chroms$Chr[i]] <- e[allbins$chr==chroms$Chr[i]] <- chroms$before[i]
    allbins$s <- s+allbins$start
    allbins$e <- e+allbins$end
    return(allbins)
}

pileup_copies <- function(regions, chroms=NULL) {
    chrs <- unique(as.character(regions$Chromosome))
    if (length(chrs)>1) { # if more than one chrom, run them separately.
        pile=NULL
        for (i in 1:length(chrs)) pile <- rbind(pile,pileup_copies(regions[as.character(regions$Chromosome)==chrs[i],],chroms))
        return(pile[order(pile$Chromosome),])
    } else {
        c=as.character(regions$Chromosome[1])
        starts=data.frame(Chromosome=c,Start=regions$Start,End=regions$Start,type='start', Cn=regions$Cn, mCn=regions$mCn)
        ends=data.frame(Chromosome=c,Start=regions$End,End=regions$End,type='end', Cn=regions$Cn, mCn=regions$mCn)
        pile <- rbind(starts,ends)
        pile <- pile[order(pile$Start),]
        cn <- c('cn0','cn1','cn2','cn3','cn4','cn5')
        try(pile$Cn[pile$Cn>5] <- 5, silent=T)
        pile$mCn[is.na(pile$mCn)] <- -1
        pile$cn0hom <- 0
        pile$cn1hom <- 0
        pile$cn2 <- 0   
        pile$cn3 <- 0
        pile$cn4 <- 0   
        pile$cn5 <- 0
        #pile$cn6 <- 0   
        pile$cn2hom <- 0   
        pile$cn3hom <- 0
        pile$cn4hom <- 0   
        pile$cn5hom <- 0
        #pile$cn6hom <- 0   
        
        n <- nrow(pile)
        thisCn <- paste('cn',pile$Cn[1],sep='')
        if (pile$mCn[1]==0) thisCn <- paste(thisCn,'hom',sep='')
        pile[[thisCn]][1:n] <- 1
        #pile$loh[1:n] <- ifelse(pile$mCn[1]==0,1,0)
        
        for (i in 2:n) {
            thisCn <- paste('cn',pile$Cn[i],sep='')
            if (pile$mCn[i]==0) thisCn <- paste(thisCn,'hom',sep='')
            pile[[thisCn]][i:n] <- pile[[thisCn]][i-1] + ifelse(pile$type[i]=='start',1,-1)            
            #if (pile$mCn[i]==0) pile$loh[i:n] <- pile$loh[i-1] + ifelse(pile$type[i]=='start',1,-1)
            pile$End[i-1] <- pile$Start[i]
        }
        pile <- pile[(pile$End > pile$Start),]
        if (!is.null(chroms)) {
            e <- 0 -> s
            for (i in 1:24) s[pile$Chromosome==as.character(chroms$Chr[i])] <- e[pile$Chromosome==as.character(chroms$Chr[i])] <- chroms$before[i]
            pile$s <- s+pile$Start
            pile$e <- e+pile$End
        }
        return(pile)
    }
}

addGenes <- function(data,genes) {
    if (nrow(data)==0) return(data)
    data$genes=''
    for (i in 1:nrow(data)) {
        ix=which(genes$chrom==as.character(data$Chromosome[i]) & genes$gtxStart<data$End[i] & genes$gtxEnd>data$Start[i])
        if (length(ix)>0) data$genes[i]=paste(sort(unique(genes$name2[ix])),collapse=', ')
    }
    return(data)
}

### Markus 2013.
### Function for summarizing alteration frquencies in X samples
TAPS_freq <- function(samples='all', outdir='frequencies', hg19=T) {
    #suppressPackageStartupMessages(library(xlsx))    
    #browser()
    sampleData <- load.txt('SampleData.csv')
    olddir <- getwd()
    if (!is.na(outdir)) {
        try(dir.create(outdir), silent=T)
        setwd(outdir)
    }
    
    subs=as.character(sampleData$Sample)
    
    if (is.null(subs)) {
        subs=thisSubdir()
        setwd('..')
    }
    
    if (samples[1]=='all') samples=rep(T,length(subs))
    if (is.logical(samples)) samples=which(samples)
    subs=subs[samples]
    nSamples=length(subs)
    
    ## Load and parse all samples
    samples <- meanCns <- NULL
    
    for (i in 1:nSamples) try( {
        #cat(subs[i],'\n')
        table <- load.txt(paste(olddir,'/',subs[i],'/',subs[i],'_segmentCN.txt',sep=''))
        table$n <- i
        table$name <- as.character(subs[i])    
        ix <- as.numeric(deChrom_ucsc(table$Chromosome)) <= 22
        table=table[!is.na(ix),]
        ix <- as.numeric(deChrom_ucsc(table$Chromosome)) <= 22
        table$meanCn <- meanCns[i] <- round(weightedMean(table$Cn_[ix], table$lengthMB[ix]),2)    
        samples <- rbind(samples,table) 
    }, silent=T)

    samples=samples[!is.na(samples$Cn),]
    regs=samples[samples$Chromosome!='chrY',]
    
    if (hg19) genes <- knownGene else genes <- knownGene_hg18
    genes$name2=genes$gAlias
    genes$chrom=genes$chr
    
    
    n=length(unique(regs$name))
    
    gaincolor='#007602'
    losscolor='#FF0000'
    lohcolor='#0015FF'
    relhighcolor='#00CE03'
    lowcolor='#FF1FE1'
    lowlohcolor='#0084FF'
    highgaincolor='#003300'
    homlosscolor='#660066'
    focalgaincolor='#333300'
    focallosscolor='#9933FF'
    
    ## Gain:
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn>2,],n,comparison='Absolute Gain (>2 copies)',color=gaincolor)
    ## Loss:
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn<2,],n,comparison='Absolute Loss (<2 copies)',color=losscolor)
    ## LOH
    aix=!is.na(regs$mCn)
    sum_regionSet(chroms, chromData, genes,regs[aix,][regs[aix,]$mCn==0,],n,comparison='LOH (minor=0)',color=lohcolor) 
    ## Relative Gain (more than 25% above the average)
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn>regs$meanCn*1.25,],n,comparison='Relative Gain (>1.25*average)',color=relhighcolor)
    ## Relative Loss (up to 3/4 of average)
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn<0.75*regs$meanCn,],n,comparison='Relative Loss (<0.75*average)',color=lowcolor)
    ## Low and LOH
    sum_regionSet(chroms, chromData, genes,regs[aix,][regs[aix,]$mCn==0 & regs[aix,]$Cn<0.75*regs$meanCn[aix],],n,comparison='Low&LOH',color=lowlohcolor)
    ## High Gain (more than 3 times the average)
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn>regs$meanCn*3,],n,comparison='High Gain (>3*average)',color=highgaincolor)
    ## Homozygous Loss
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn==0,],n,comparison='Homozygous Loss (0 copies)',color=homlosscolor)
    ## Focal Gain
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn>regs$meanCn*2 & regs$lengthMB<1,],n,comparison='Focal Gain (<1Mb,>2*average)',color=focalgaincolor)
    ## Focal Loss
    sum_regionSet(chroms, chromData, genes,regs[regs$Cn<regs$meanCn*.67 & regs$lengthMB<1,],n,comparison='Focal Loss (<1Mb,<0.67*average)',color=focallosscolor)
    setwd(olddir)
}

### Generalized comparing function used by the above function 
sum_regionSet <- function(chroms, chromData, genes, 
                          regs, n, 
                          comparison='', 
                          color='#000000') {    
    
    p_cutoff <- 0.05
    freq_cutoff <- 20
    
    # Generate frequency plots            
    pile <- pileup(regs); pile$Percent <- 100*pile$Count / n
    save.txt(addGenes(pile,genes),file=paste(comparison,'frequencies','txt',sep='.'))
    
    ### Total frequency plot    
    #quartz(file=paste(comparison,'freq_dif',name1,name2,'png',sep='.'),width=15,height=4.5,dpi=300,type='png')
    png(filename=paste(comparison,'png',sep='.'),width=3000,height=800)
    ylim <- c( -15, 110)
    plot(1,1,type='n',
         bty='n', ann=T, mar=c(0,0,0,0), oma=c(0,0,0,0),
         main=paste(comparison,'frequency'), cex.main=2,
         xlab=NA,#'',#Genomic position (??????=significance)',
         ylab=NA,#'Alteration frequency difference (% units)', cex.lab=2,
         xlim = c(0,sum(as.numeric(chroms$length))),
         ylim = ylim,
         yaxt="n",
         xaxt="n"
         )
    #mtext(text="Allelic imbalance",side=2,line=4,las=3,cex=2)
    #mtext(text="Coverage, all segments",side=1,line=4,cex=2,las=1)
    
    if (nrow(pile)>0) {
        data <- pile
        data$e <- 0 -> data$s
        for (i in 1:23) data$s[data$Chromosome==as.character(chroms$Chr[i])] <- chroms$before[i] -> data$e[data$Chromosome==as.character(chroms$Chr[i])]
        data$s <- data$s+data$Start
        data$e <- data$e+data$End
        rect( 
            xleft=data$s, xright=data$e,
            ybottom=0,ytop=data$Percent,
            col=color,   
            border=NA,
        )
    }
    ## The chrom labels
    text(
        x=chroms$before+chroms$length/2,
        y= -15,
        cex=2,
        label=substr(chroms$chr,4,6)
    )
    
    ## The chromosomes
    rect(xleft=chromData$genomeStart,xright=chromData$genomeEnd,
         ybottom= -5-chromData$thickness,ytop= -5+chromData$thickness,
         col=paste(chromData$col,'80',sep=''),
         lty=0
    )
    
    ## The vertical lines
    segments(
        x0=c(chroms$before,sum(as.numeric(chroms$length))),x1=c(chroms$before,sum(as.numeric(chroms$length))),
        y0= -15,y1=100,                
        col='#000000',
        lwd=1
    )
    
    # Horizontal lines
    segments(
        x0=0,x1=sum(as.numeric(chroms$length)),
        y0=seq(0,80,10),y1=seq(0,80,10),                
        col='#00000070',
        lwd=1,lty=3
    )
    
    ## Y axis ticks
    axis(2,at=seq(0,100,20), 
         labels=abs(seq(0,100,20)), pos=0,
         cex.axis=2,tck=-0.01, las=2)
    axis(4,at=seq(0,100,20), 
         labels=abs(seq(0,100,20)), pos=sum(as.numeric(chroms$length)),
         cex.axis=2,tck=-0.01, las=2)
    dev.off()  
}

########### Group comparisons
TAPS_compare <- function(grp1, grp2, name1='1', name2='2', outdir='frequencies_comp', hg19=T) {

    #suppressPackageStartupMessages(library(xlsx))    

    sampleData=load.txt('SampleData.csv')
    
    subs=as.character(sampleData$Sample)
    
#     if (is.null(subs)) {
#         subs=thisSubdir()
#         subs=subs[subs!='frequencies' & subs!='frequencies_comp']
#         setwd('..')
#     }
    
    olddir <- getwd()
    if (!is.na(outdir)) {
        try(dir.create(outdir), silent=T)
        setwd(outdir)
    }
    

    
    if (is.logical(grp1)) grp1=which(grp1)
    if (is.logical(grp2)) grp2=which(grp2)
    if (is.numeric(grp1)) grp1=subs[grp1]
    if (is.numeric(grp2)) grp2=subs[grp2]
    grps=unique(c(grp1,grp2))
    
    ## Load and parse all samples
    samples <- meanCns <- NULL
    nSamples=length(grps)
    
    for (i in 1:nSamples) {
        table <- load.txt(paste(olddir,'/',grps[i],'/',grps[i],'_segmentCN.txt',sep=''))
        table$n <- i
        table$name <- as.character(grps[i])    
        ix <- as.numeric(deChrom_ucsc(table$Chromosome)) <= 22
        table$meanCn <- meanCns[i] <- round(weightedMean(table$Cn_[ix], table$lengthMB[ix]),2)    
        samples <- rbind(samples,table) 
    }
    
    samples=samples[!is.na(samples$Cn),]
    samples=samples[samples$Chromosome!='chrY',]
    
    if (hg19) genes <- knownGene else genes <- knownGene_hg18
    genes$name2=genes$gAlias
    genes$chrom=genes$chr
    
    regs1=samples[samples$name %in% grp1,]
    regs2=samples[samples$name %in% grp2,]
    n1=length(grp1)
    n2=length(grp2)
    
    gaincolor='#007602'
    losscolor='#FF0000'
    lohcolor='#0015FF'
    relhighcolor='#00CE03'
    lowcolor='#FF1FE1'
    lowlohcolor='#0084FF'
    highgaincolor='#003300'
    homlosscolor='#660066'
    focalgaincolor='#333300'
    focallosscolor='#9933FF'
    
    ## Gain:
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn>2,],regs2[regs2$Cn>2,],n1,n2,comparison='Absolute Gain (>2 copies)',name1,name2,color=gaincolor)
    ## Loss:
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn<2,],regs2[regs2$Cn<2,],n1,n2,comparison='Absolute Loss (<2 copies)',name1,name2,color=losscolor)
    ## LOH
    aix1=!is.na(regs1$mCn); aix2=!is.na(regs2$mCn)
    compare_regionSet(chroms, chromData, genes,regs1[aix1,][regs1[aix1,]$mCn==0,],regs2[aix2,][regs2[aix2,]$mCn==0,],n1,n2,comparison='LOH (minor=0)',name1,name2,color=lohcolor) 
    ## Relative Gain (more than 25% above the average)
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn>regs1$meanCn*1.25,],regs2[regs2$Cn>regs2$meanCn*1.25,],n1,n2,comparison='Relative Gain (>1.25*average)',name1,name2,color=relhighcolor)
    ## Relative Loss (up to 3/4 of average)
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn<0.75*regs1$meanCn,],regs2[regs2$Cn<0.75*regs2$meanCn,],n1,n2,comparison='Relative Loss (<0.75*average)',name1,name2,color=lowcolor)
    ## Low and LOH
    compare_regionSet(chroms, chromData, genes,regs1[aix1,][regs1[aix1,]$mCn==0 & regs1[aix1,]$Cn<0.75*regs1$meanCn[aix1],],regs2[aix2,][regs2[aix2,]$mCn==0 & regs2[aix2,]$Cn<0.75*regs2$meanCn[aix2],],n1,n2,comparison='Low&LOH',name1,name2,color=lowlohcolor)
    ## High Gain (more than 3 copies above the average)
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn>regs1$meanCn*3,],regs2[regs2$Cn>regs2$meanCn*3,],n1,n2,comparison='High Gain (>3*average)',name1,name2,color=highgaincolor)
    ## Homozygous Loss
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn==0,],regs2[regs2$Cn==0,],n1,n2,comparison='Homozygous Loss (0 copies)',name1,name2,color=homlosscolor)
    ## Focal Gain
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn>regs1$meanCn*2 & regs1$lengthMB<1,],regs2[regs2$Cn>regs2$meanCn*2 & regs2$lengthMB<1,],n1,n2,comparison='Focal Gain (<1Mb,>2*average)',name1,name2,color=focalgaincolor)
    ## Focal Loss
    compare_regionSet(chroms, chromData, genes,regs1[regs1$Cn<regs1$meanCn*.67 & regs1$lengthMB<1,],regs2[regs2$Cn<regs2$meanCn*.67 & regs2$lengthMB<1,],n1,n2,comparison='Focal Loss (<1Mb,<0.67*average)',name1,name2,color=focallosscolor)
    setwd(olddir)
}

### Generalized comparing function: 
compare_regionSet <- function(chroms, chromData, genes, 
                              regs1, regs2, 
                              n1, n2, 
                              comparison='', 
                              name1='1', name2='2', 
                              color='#000000') {    
    p_cutoff <- 0.05
    freq_cutoff <- 0
    
    difcolor <- color
    color <- colorRampPalette(c("#FFFFFF",color),space="rgb")(4)[2]
    
    # Generate frequency plots            
    pile1 <- pileup(regs1); pile1$Percent <- 100*pile1$Count / n1
    pile2 <- pileup(regs2); pile2$Percent <- 100*pile2$Count / n2
    dif <- pvals(pileup_dif(regs1,regs2),n1,n2)
    
    save.txt(addGenes(dif,genes),file=paste(comparison,name1,name2,'txt',sep='.'))
    save.txt(addGenes(dif[dif$p.value<p_cutoff & abs(dif$percent)>freq_cutoff,],genes),file=paste('significant',comparison,name1,name2,'txt',sep='.'))
    
    
    ### Total and Differential frequency plot    
    #quartz(file=paste(comparison,'freq_dif',name1,name2,'png',sep='.'),width=15,height=4.5,dpi=300,type='png')
    png(filename=paste(comparison,'allFreqs',name1,name2,'png',sep='.'),width=3000,height=800)
    ylim <- c(-100, 110)
    plot(1,1,type='n',
         bty='n', ann=T, mar=c(0,0,0,0), oma=c(0,0,0,0),
         main=paste(comparison,' frequency difference, ', name1, ' (', n1, ') vs ', name2, ' (', n2, ')', sep=''), cex.main=2,
         xlab=NA,#'',#Genomic position (??????=significance)',
         ylab=NA,#'Alteration frequency difference (% units)', cex.lab=2,
         xlim = c(0,sum(as.numeric(chroms$length))),
         ylim = ylim,
         yaxt="n",
         xaxt="n"
         )
    #mtext(text="Allelic imbalance",side=2,line=4,las=3,cex=2)
    #mtext(text="Coverage, all segments",side=1,line=4,cex=2,las=1)
    mhcolor=colorRampPalette(c("#FFFFFF",color),space="rgb")(5)[3]
    
    if (nrow(pile1)>0) {
        ## Group 1
        data <- pile1
        data$e <- 0 -> data$s
        for (i in 1:23) data$s[data$Chromosome==as.character(chroms$Chr[i])] <- chroms$before[i] -> data$e[data$Chromosome==as.character(chroms$Chr[i])]
        data$s <- data$s+data$Start
        data$e <- data$e+data$End
        rect( 
            xleft=data$s, xright=data$e,
            ybottom=0,ytop=data$Percent,
            col=color,   
            border=NA,
        )
    }
    if (nrow(pile2)>0) {
        ## Group 2
        data <- pile2
        data$e <- 0 -> data$s
        for (i in 1:23) data$s[data$Chromosome==as.character(chroms$Chr[i])] <- chroms$before[i] -> data$e[data$Chromosome==as.character(chroms$Chr[i])]
        data$s <- data$s+data$Start
        data$e <- data$e+data$End
        rect( 
            xleft=data$s, xright=data$e,
            ybottom=0,ytop= -data$Percent,
            col=color,   
            border=NA,
        )
    }
    if (nrow(dif)>0) {
        ### The difference
        data <- dif
        data$e <- 0 -> data$s
        for (i in 1:23) data$s[data$Chromosome==as.character(chroms$Chr[i])] <- chroms$before[i] -> data$e[data$Chromosome==as.character(chroms$Chr[i])]
        data$s <- data$s+data$Start
        data$e <- data$e+data$End
        data <- data[order(data$percent,decreasing=T),]
        
        #difcolor <- colorRampPalette(c("#FFFFFF",color),space="rgb")(4)[2]
        
        rect( 
            xleft=data$s, xright=data$e,
            ybottom=0,ytop=data$percent,
            #ybottom=data$conf_low,ytop=data$conf_high,
            col=difcolor,   
            border=NA
        )
        
        ## The significant
        data=data[data$p.value<p_cutoff & abs(data$percent)>freq_cutoff,]
        if (nrow(data)>0) {
            rect( 
                xleft=data$s, xright=data$e,
                ybottom= -100,ytop=-95,
                col='#000000',   
                border=NA
            )
        }
        ## The OR
        #for (i in 1:23) {
        #    temp=data[data$Chr==as.character(chroms$chr[i]),]
        #    if (nrow(temp)>5) {
        #        temp$odds.ratio[temp$odds.ratio<1 & temp$odds.ratio>0]=1/temp$odds.ratio[temp$odds.ratio<1 & temp$odds.ratio>0]
        #        temp=temp[order(temp$odds.ratio,decreasing=T),][1,]
        #        text(x=chroms$before[i]+mean(temp$Start,temp$End), y=-105, cex=1.3, label=round(temp$odds.ratio,1))                
        #    }
        #}
    }
    axis(2,at=c(-97),#,-105), 
         labels=c('sign.'),#, 'BestOR'), 
         pos=0,
         cex.axis=1.5,tck=0,las=2)
    
    ## The chrom labels
    text(
        x=chroms$before+chroms$length/2,
        y= 105,
        cex=2,
        label=substr(chroms$chr,4,6)
    )
    
    ## The chromosomes
    rect(xleft=chromData$genomeStart,xright=chromData$genomeEnd,
         ybottom=95-chromData$thickness,ytop=95+chromData$thickness,
         col=paste(chromData$col,'80',sep=''),
         lty=0
    )
    
    ## The vertical lines
    segments(
        x0=c(chroms$before,sum(chroms$length)),x1=c(chroms$before,sum(chroms$length)),
        y0=-100,y1=100,                
        col='#000000',
        lwd=1
    )
    
    # Horizontal lines
    segments(
        x0=0,x1=sum(chroms$length),
        y0=seq(-80,80,10),y1=seq(-80,80,10),                
        col='#00000070',
        lwd=1,lty=3
    )
    segments(
        x0=0,x1=sum(chroms$length),
        y0=0,y1=0,                
        col='#00000090',
        lwd=1,lty=1
    )
    
    ## Y axis ticks
    axis(2,at=seq(-80,80,20), 
         labels=abs(seq(-80,80,20)), pos=0,
         cex.axis=2,tck=-0.01, las=2)
    axis(4,at=seq(-80,80,20), 
         labels=abs(seq(-80,80,20)), pos=sum(chroms$length),
         cex.axis=2,tck=-0.01, las=2)
    
    ## Group names
    mtext(c(paste(name1,'(%)'), '(dif)', paste(name2,'(%)')),2,at=c(65,0,-65),cex=2, col=c(color,mhcolor,color))
    
    dev.off()      
}



#------------------------------------------------------------
# Author: Sebastian DiLorenzo , 2013
#------------------------------------------------------------

#------------------------------------------------------------
# Added TAPS plot functionality
#------------------------------------------------------------

# #######  ##     ## ######## ########  ##     ## #### ######## ##      ## ########  ##        #######  ######## 
###     ## ##     ## ##       ##     ## ##     ##  ##  ##       ##  ##  ## ##     ## ##       ##     ##    ##    
###     ## ##     ## ##       ##     ## ##     ##  ##  ##       ##  ##  ## ##     ## ##       ##     ##    ##    
###     ## ##     ## ######   ########  ##     ##  ##  ######   ##  ##  ## ########  ##       ##     ##    ##    
###     ##  ##   ##  ##       ##   ##    ##   ##   ##  ##       ##  ##  ## ##        ##       ##     ##    ##    
###     ##   ## ##   ##       ##    ##    ## ##    ##  ##       ##  ##  ## ##        ##       ##     ##    ##    
# #######     ###    ######## ##     ##    ###    #### ########  ###  ###  ##        ########  #######     ##                                                            

OverviewPlot <- function(chr,start,end,int,ai,hg18,mchr,mpos,mval,schr,spos,sval,name='',xlim=c(-1,1.5),ylim=c(0,1),MAPD,MHOF)
{
    if(hg18==T)
        {
        chroms=chroms_hg18
        chromData=chromData_hg18
        chroms$length=as.numeric(chroms$length)
        }

    #ideogram=getIdeogram()
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
    
    size=rep(0.2,length(chr)); size[int<xlim[1]]=0; size[int>xlim[2]]=0
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
    
    split.screen(as.matrix(data.frame(left=c(rep(0.015,3)),
                                      right=c(rep(0.987,3)),
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
        this <- chroms[chroms$c==c,]
        #Extract that chromosomes information
        ix <- chr==as.character(this$chr)
        x <- chr=='chrX'
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
        plot(c(int[!(ix|x)],int[ix]),c(ai[!(ix|x)],ai[ix]),
             pch=c(pch[!(ix|x)],pch[ix]),
             cex=c(size[!(ix|x)],size[ix]),
             main = "",
             xlab = '',
             ylab = "",
             col = c(col[!(ix|x)],col[ix]),
             xlim=xlim,ylim=c(0,1),
             axes=F)
        
        #grid lines
        segments(
            x0=c(-0.5,0,0.5),x1=c(-0.5,0,0.5),
            y0=-0.05,y1=1.035,                
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
        if(c %in% c(seq(from=17,to=24,by=1))) axis(side=1,cex.axis=0.6,at=seq(from=-0.5,to=1,by=0.5))
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
            mtext("Average log-ratio",side=1,line=1.5)
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
    par(lend=1)
    
    #Remove chrY,M,XY
    mval[mchr=='chrXY'|mchr=='chrM']=NA    
    # mval[mchr=='chrY'|mchr=='chrXY'|mchr=='chrM']=NA    
    
    #Calculate previous distance of whole genome so the next chromsome is added
    #at the correct coordinate
    pre=rep(NA,24)
    for (c in 1:24)
    {
        this <- chroms[chroms$c==c,]
        if(nrow(this) == 0) next
        ix <- chr==as.character(this$chr)
        mix <- mchr==as.character(this$chr) #& mval>=(-2)
        six <- schr==as.character(this$chr)
        
        #Predefine ymin ymax and sequence between them
        ymin=floor(2*min(int,na.rm=T))/2
        if(ymin > -1)
        {
            ymin = -1
        }
        ymax=ceiling(2*max(int,na.rm=T)+1)/2
        if(ymax < 1)
        {
            ymax = 1
        }
        seqminmax=seq(-2,2,by=0.5)
        
        if(c>1)
        {
            prelength=prelength+chroms$length[chroms$c==(c-1)]
        } else
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
         col = '#00000003',
         xaxt="n",
         axes=F,
         ylim = c(-2,2),
         xlim = c(0,max(mpos))
    )
    
    #Add axis to the left & right of signal
    axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
    axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=sum(chroms$length),las=1)
    mtext("Log-ratio",side=2,line=0)
    
    #Add grey segments
    segments(
        y0=c(-0.5,0,0.5),y1=c(-0.5,0,0.5),
        x0=0,x1=max(mpos),                
        col='#D3D3D360',
        lwd=1)
    
    #Add colored segment information for each chromosome, red to blue gradient.
    for(c in 1:24)
    {
        this <- chroms[chroms$c==c,]
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
        x0=c(chroms$before,sum(chroms$length)),x1=c(chroms$before,sum(chroms$length)),
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
    
    #Remove M,XY
    sval[mchr=='chrXY'|mchr=='chrM']=NA    
    
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
    axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X",'Y'),cex.axis=0.55,lty=0)#,tck=0,col.ticks='#00000000')
    axis(side=4,tck=-0.04,at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=sum(chroms$length),las=1) #
    mtext("Allele frequency",side=2,line=0)
    mtext("Chromosomes",side=1,line=1.5,adj=0.4)
    
    #Add a bar between chromosomes to distinguish them
    segments(
        x0=c(chroms$before,sum(chroms$length)),x1=c(chroms$before,sum(chroms$length)),
        y0=0,y1=100,                
        col='#000000',
        lwd=1)
    
    
    #Close all the opened split.screens and release the figure
    close.screen(all.screens=T)
    dev.off()
    
    file.copy(paste(name,'_overview.jpg',sep=''),paste("../",name,'_overview.jpg',sep=''),overwrite=T)
    
}


###    ##    ###    ########  ##    ##  #######  ######## ##    ## ########  ######## 
###   ##    ## ##   ##     ##  ##  ##  ##     ##    ##     ##  ##  ##     ## ##       
###  ##    ##   ##  ##     ##   ####   ##     ##    ##      ####   ##     ## ##       
######    ##     ## ########     ##    ##     ##    ##       ##    ########  ######   
###  ##   ######### ##   ##      ##    ##     ##    ##       ##    ##        ##       
###   ##  ##     ## ##    ##     ##    ##     ##    ##       ##    ##        ##       
###    ## ##     ## ##     ##    ##     #######     ##       ##    ##        ######## 


# ######  ##     ## ########   #######  ##     ##  ######  
###    ## ##     ## ##     ## ##     ## ###   ### ##    ## 
###       ##     ## ##     ## ##     ## #### #### ##       
###       ######### ########  ##     ## ## ### ##  ######  
###       ##     ## ##   ##   ##     ## ##     ##       ## 
###    ## ##     ## ##    ##  ##     ## ##     ## ##    ## 
# ######  ##     ## ##     ##  #######  ##     ##  ######  

#Function for generating individual plots for TAPS_call. Gives whole genome, log-ratio, cytoband and allele frequency.
karyotype_chroms <- function(chr,start,end,int,ai,hg18,mchr,mpos,mval,schr,spos,sval,name='',xlim=c(-2,2),ylim=0:1,MAPD,MHOF)  
{   
    if(hg18==T)
    {
    chroms=chroms_hg18
    chromData=chromData_hg18
    }
    #Get ideogram  
    #ideogram=getIdeogram()
    
    #Set color gradient for p and q arm of chromosome
    colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
    colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")
    
    #filter the data (remove NA) and set some standard parameters
    #ai[is.na(ai)]=0
    #aix=ai!=0
    #chr=chr[aix]
    #start=start[aix]
    #end=end[aix]
    #int=int[aix]
    #ai=ai[aix]
    pos <- (start+end)/2
    length=end-start
    
    #Set sizes for circles in whole genome plot
    size=rep(0.5,length(chr))#; size[chr=='chrX']=0.2
    #size[length>2000000]=2
    #size[length>5000000]=3
    #size[length>10000000]=4
    
    #Loop over each chromosome. (Not chromosomeY as of writing)     
    for (c in 1:24) 
    {
        
        #Select the chromosome
        this <- chroms[chroms$c==c,]
        #Skip current chromosome if it doesn't exist.
        if(!any(as.character(this$chr) == chr)) next()
        if(nrow(this) == 0) next
        
        #Initialize jpeg
        jpeg(paste(name,'_karyotype.',this$chr,'.jpg',sep=''),width=11.7,height=8.3,units="in",res=300)
        
        #split plot into desired formation
        split.screen(as.matrix(data.frame(left=c(rep(0.05,4)),
                                          right=c(rep(1,4)),
                                          bottom=c(0.05,0.20,0.25,0.50),
                                          top=c(0.20,0.25,0.40,0.95)))) #screen 1,2,3
        
        #Screen configuration overview
        #------------------------------------------------------------
        #|                                                          |
        #|                                                          |
        #|                                                          |
        #|                                                          |
        #|                            4                             |
        #|                                                          |
        #|                                                          |
        #|                                                          |
        #|                                                          |
        #------------------------------------------------------------
        #|                            3                             |
        #------------------------------------------------------------
        #|                            2                             |
        #------------------------------------------------------------
        #|                            1                             |
        #------------------------------------------------------------ 
        
        #------------------------------------------------------------
        #Top part of the plot
        #------------------------------------------------------------
        
        #index of the selected chromsome to the chr object
        ix <- chr==as.character(this$chr)
        
        notix1 = chr=="chrX"
        notix = !ix + notix1
        
        #Create an index of colors relating to the positions and lengths on this chromosome
        col <- rep('#46464620',length(chr))       
        col[ix & (pos < this$mid)] <- paste(colors_p(sum(ix & (pos < this$mid))), '70', sep='')  
        col[ix & (pos > this$mid)] <- paste(colors_q(sum(ix & (pos > this$mid))), '70', sep='')
        
        #Go to screen 4
        screen(4)

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
                 xlim = xlim,
                 ylim = ylim)
        } else
        {
            plot(int[notix],ai[notix],
                 pch=16,
                 cex=size[notix],
                 main = "",
                 xlab = "",
                 ylab = "",
                 axes=F, #Remove axis
                 col = col[notix],
                 xlim = xlim,
                 ylim = ylim)
        }
        
        xix=chr=="chrX"
        if(c!=23)
        {
            points(int[xix],ai[xix],pch=4,cex=0.5,col='#46464620')
        } else
        {
            points(int[xix],ai[xix],pch=4,cex=0.5,col=col[ix]) 
        }
        
        #Insert Y and X axis
        axis(side=2,cex.axis=0.6,tck=0.925,col.ticks='#80808040',pos=xlim[1],at=seq(from=0,to=1,by=0.1),las=1)
        axis(side=1,cex.axis=0.6,tck=0.925,col.ticks='#80808040',pos=0,at=seq(from=xlim[1],to=xlim[2],by=0.1))
        
        #Titles,date/time and axis labels
        mtext(text="Average log-ratio",side=1,line=1,cex=1)
        mtext(text="Allelic imbalance",side=2,line=0,cex=1)
        mtext(paste(name,"\nChromsome ",c,sep=""),side=3)
        mtext(paste('MAPD: ',MAPD,', MHOF: ',MHOF,'\n',format(Sys.time(),"%Y-%m-%d %H:%M")),side=3,cex=0.6,adj=0.95)
        
        #------------------------------------------------------------
        #Middle - Signal
        #------------------------------------------------------------
        
        #Go to screen 3
        screen(3)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        par(lend=1)
        
        #Select the correct chromosome and remove stuff lower than -1 
        mix <- mchr==as.character(this$chr) #& mval>(-1)
        
        #Predefine ymin ymax and sequence between them
        ymin=floor(2*min(int,na.rm=T))/2
        if(ymin > -1)
        {
            ymin = -1
        }
        ymax=ceiling(2*max(int,na.rm=T)+1)/2
        if(ymax < 1)
        {
            ymax = 1
        }
        #ymin=floor(min(int[ix]))-0.5
        #if(ymin > -1)
        #{
        #    ymin = -1
        #}
        #ymax=ceiling(max(int[ix]))+0.5
        #if(ymax < 1)
        #{
        #    ymax = 1
        #}
        seqminmax=seq(-2,2,by=0.5)
        
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
             ylim = c(-2,2))
        
        #Add colored segments based on the log-ratio data
        segments(x0=start[ix],x1=end[ix],
                 y0=int[ix],y1=int[ix],                
                 col=col,
                 lwd=2)
        
        #Add X and Y axis. The (side=4 and =1) axis is just a black line showing where the data ends
        # 2:left
        axis(side=2,tck=0.926,col.ticks='#80808040',at=seqminmax, #seq(from=-1,to=1.5,by=0.5),
             #labels=c("-1","-.5","0",".5","1","1.5"),
             cex.axis=0.6,pos=0,las=1)
        # 4: right
        axis(side=4,tck=0,pos=this$length,at=seqminmax,
             #labels=c("-1","-.5","0",".5","1","1.5"),
             cex.axis=0.6,las=1)
        # 1: below
        # axis(side=1,tck=1.041,col='#80808040',at=seq(5e6,this$length,by=5e6),labels=FALSE,cex.axis=0.6,pos=ymin)
        segments(x0=seq(5e6,this$length,by=5e6),x1=seq(5e6,this$length,by=5e6),y0=-2,y1=2,col='#80808040')
        #axis(side=1,tck=0,pos=ymin,
        #     cex.axis=0.6,las=1)
        # 3: top
        #axis(side=3,tck=0,pos=ymax,
        #     cex.axis=0.6,las=1)
        
        #Add Y axis label
        mtext("Log-ratio",side=2,line=-0.8)
        
        #------------------------------------------------------------
        #Cytobands
        #------------------------------------------------------------
        
        #Go to screen 2
        screen(2)
        
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(0.5,0.25,0))
        
        #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
        plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,0.3),xlim=c(0,max(mpos[mix])))
        
        #Add Y axis label
        mtext(text=this$chr,side=2,las=1,line=-1.2)
        
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
            if(T) #(chromData$chromEnd[dix][j] - chromData$chromStart[dix][j]) > 4*10^6)
            {
                # Text is added at x position (all of previous chromosome) + (middle of region)
                # The "srt" flips the text 90 degrees
                text(x=(chromData$chromStart[dix][j]+((chromData$chromEnd[dix][j]-chromData$chromStart[dix][j])/2)),
                     y=0.15,labels=chromData$name[dix][j], srt=90,cex = 0.3,xpd=T)
            }
        }
        
        #------------------------------------------------------------
        #Bottom - AI
        #------------------------------------------------------------ 
        
        #Go to screen 1
        screen(1)
        
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
        
        #Add X and Y axis as well as a line to the rightmost of the data (side=4 and =1)
        axis(side=2,tck=0.926,col.ticks='#80808080',at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=0,las=1)
        axis(side=4,tck=0,at=seq(from=0,to=1,by=0.2),cex.axis=0.6,pos=this$length,las=1)
        axis(side=1,tck=0.925,col='#80808000',col.ticks='#80808040',at=seq(5e6,this$length,by=5e6),
             labels=paste(seq(5,floor(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
        #axis(side=1,tck=0.926,col='#80808040',at=seq(5e6,this$length,by=5e6),
        #     labels=FALSE,cex.axis=0.6,pos=ymin)
#         axis(side=1,tck=0,pos=0,at=c(0,max(mpos)),
#              cex.axis=0.6,las=1)
#         axis(side=3,tck=0,pos=1,at=c(0,max(mpos)),
#              cex.axis=0.6,las=1)
#         
        #Add X and Y label
        mtext("Allele frequency",side=2,line=-0.8)
        mtext("Position (Mb)",side=1,line=1)
        
        #Close all the opened split.screens and release the figure
        close.screen(all.screens=T)
        dev.off()
    }
}

###    ##    ###    ########  ##    ##  #######  ######## ##    ## ########  ######## 
###   ##    ## ##   ##     ##  ##  ##  ##     ##    ##     ##  ##  ##     ## ##       
###  ##    ##   ##  ##     ##   ####   ##     ##    ##      ####   ##     ## ##       
######    ##     ## ########     ##    ##     ##    ##       ##    ########  ######   
###  ##   ######### ##   ##      ##    ##     ##    ##       ##    ##        ##       
###   ##  ##     ## ##    ##     ##    ##     ##    ##       ##    ##        ##       
###    ## ##     ## ##     ##    ##     #######     ##       ##    ##        ######## 


# ######  ##     ## ########   #######  ##     ##  ######      ######  ##    ## 
###    ## ##     ## ##     ## ##     ## ###   ### ##    ##    ##    ## ###   ## 
###       ##     ## ##     ## ##     ## #### #### ##          ##       ####  ## 
###       ######### ########  ##     ## ## ### ##  ######     ##       ## ## ## 
###       ##     ## ##   ##   ##     ## ##     ##       ##    ##       ##  #### 
###    ## ##     ## ##    ##  ##     ## ##     ## ##    ##    ##    ## ##   ### 
# ######  ##     ## ##     ##  #######  ##     ##  ######      ######  ##    ## 

karyotype_chromsCN <- function(chr,start,end,int,ai,Cn,mCn,hg18,mchr,mpos,mval,schr,spos,sval,t,name='',xlim=c(-1.02,1.82),ylim=0:1, maxCn=8,parameters)  
{   
    if(hg18==T) {
        chroms=chroms_hg18
        chromData=chromData_hg18
    }
    #Get ideogram
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
    Cn=Cn[aix]
    mCn=mCn[aix]
    pos <- (start+end)/2
    length=end-start
    labels=paste(Cn[!is.na(mCn)],mCn[!is.na(mCn)],sep='m')
    
    #Set sizes for circles in whole genome plot
    size=rep(1,length(chr))
    size[length>2000000]=2
    size[length>5000000]=3
    size[length>10000000]=4
    
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
        this <- chroms[chroms$c==c,]
        
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
        
        notix1 = chr=="chrX"
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
        
        if(c!=23) {
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
                 xlim = c(-1,1.5),
                 ylim = ylim)
        } else {
            plot(int[notix],ai[notix],
                 pch=16,
                 cex=size[notix],
                 main = "",
                 xlab = "",
                 ylab = "",
                 axes=F, #Remove axis
                 col = col[notix],
                 xlim = c(-1,1.5),
                 ylim = ylim)
        }
        
        xix=chr=="chrX"
        if(c!=23) {
            points(int[xix],ai[xix],pch=4,col='#46464680')
        } else {
            points(int[xix],ai[xix],pch=4,col=col[ix]) 
        }
        
        #Insert Y and X axis
        axis(side=2,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),las=1)
        axis(side=1,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=-1,to=1.5,by=0.1))
        
        #Add cn and mcn labels
        text(variants_data[,1],variants_data[,2],
             labels=variants,
             col='black',
             cex=0.75)
        
        #Titles and axis labels
        mtext(text="Average log-ratio",side=1,line=1,cex=1)
        mtext(text="Allelic imbalance",side=2,line=1,cex=1)
        
        #Check if name of sample is too long to be graphically pleasing and cut it if that is the case
        if(nchar(name)>12)
        {
            mtext(paste("Sample: ",substring(name,1,12),"\n Chromsome ",c,sep=""),side=3)
        }
        else
        {
            mtext(paste("Sample: ",name,"\n Chromsome ",c,sep=""),side=3)
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
        par(lend=1)
        
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
        

        #Add minor copynumber segments
        segments(x0=start[ix],x1=end[ix],
                 y0=mCn[ix],y1=mCn[ix],                
                 col='#BBBBBB',
                 lwd=5)
        
        #Add total copynumber segments
        segments(x0=start[ix],x1=end[ix],
                 y0=Cn[ix],y1=Cn[ix],                
                 col=col,
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
        par(lend=1)
        
        #Select the correct chromosome and remove stuff lower than -1 
        mix <- mchr==as.character(this$chr) #& mval>(-1)
        
        #Predefine ymin ymax and sequence between them
        ymin=floor(min(int[ix]))-0.5
        if(ymin > -1)
        {
            ymin = -1
        }
        ymax=ceiling(max(int[ix]))+0.5
        if(ymax < 1)
        {
            ymax = 1
        }
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
        mtext("log-ratio",side=2,line=0.3)
        
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
             labels=paste(seq(5,floor(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
        
        #Add X and Y label
        mtext("Allele frequency",side=2,line=0.3)
        mtext("Position (Mb)",side=1,line=1)

        #text describing the parameters used to generate the plot
        mtext(text=parameters,side=1,line=6,cex=0.7,adj=0.95)
        
        #Close all the opened split.screens and release the figure
        close.screen(all.screens=T)
        dev.off()
    }
}

#########    ###    ########   ######     ########  ########  ######   ####  #######  ##    ## 
#   ##      ## ##   ##     ## ##    ##    ##     ## ##       ##    ##   ##  ##     ## ###   ## 
#   ##     ##   ##  ##     ## ##          ##     ## ##       ##         ##  ##     ## ####  ## 
#   ##    ##     ## ########   ######     ########  ######   ##   ####  ##  ##     ## ## ## ## 
#   ##    ######### ##              ##    ##   ##   ##       ##    ##   ##  ##     ## ##  #### 
#   ##    ##     ## ##        ##    ##    ##    ##  ##       ##    ##   ##  ##     ## ##   ### 
#   ##    ##     ## ##         ######     ##     ## ########  ######   ####  #######  ##    ## 

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#TAPS_region() ; a function for zooming in on a specific region and showing the regions genes.
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


TAPS_region <- function(directory=NULL,chr,region,hg18=F)
{
    if (is.null(directory))
    {
        cat("No directory supplied, using working directory.")
        directory = "."
        #cat("You have not assigned a directory containing one or more folders of samples for TAPS_plot to execute. \n")
        #cat("Example: \"/user/mysamples/\" or, to run it in your current working directory, \".\" \n")
        #directory = readline("Please supply such a directory now: ")
    }
    
    #Set working directory and make that directorys name the samplename. If it is too long, shorten it.
    setwd(directory)
    subs <- getSubdirs()
    if (is.null(subs)) 
    {                                         ## check samples = subdirectories or a single sample = current directory
        subs=thisSubdir()
        setwd('..')
    }
    
    #store chr value in nchr because i made the newbie mistake of having a second variable also
    #named chr.
    nchr = chr  
    
    #Loop over all subdirectories, should they exist.
    for (i in 1:length(subs)) try( 
{
    setwd(subs[i])
    name <- subs[i]
    cat(' ..plotting', subs[i],"\n")
    
    if(nchar(name)>12)
    {
        name = substring(name,1,12)
    }
    
    #From TAPS_plot():
    #save(regs,Log2,alf,segments,file="TAPS_plot_output.Rdata")  
    
    #Check if TAPS_plot_output.Rdata exists. It must for the plot function to run.
    #In this file are the objects needed to plot.
    if(file.exists("TAPS_plot_output.Rdata")==F)
    {
        cat("TAPS_region could not find a TAPS_plot_output.Rdata file in the sample directory. \n
       You must run TAPS_plot before TAPS_region. \n")
    }
    else
    {
        load("TAPS_plot_output.Rdata")
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
    
    #Rename parameters of TAPS_plot_output so we can re-use other plots as template.
    chr =regs$chr
    start = regs$start  
    end = regs$end
    int = regs$logs
    ai = regs$scores
    #ideogram=NULL
    mchr = as.character(Log2$Chromosome)
    mpos =  Log2$Start
    mval =Log2$Value
    schr = as.character(alf$Chromosome) 
    spos = alf$Start
    sval =alf$Value
    xlim=c(-1,2)
    ylim=c(0,1)
    
    #check if hg18 is true
    #then load hg18 knownGene list
    if(hg18==T)
    {
        kg = knownGene_hg18
        chroms = chroms_hg18
        chromData = chromData_hg18
    }
    else
    {
        kg = knownGene
    }
    
    
    #
    #Here is where pretty normal plotting procedure starts
    #
    
    #Get ideogram  
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
    
    #Set sizes for circles in whole genome plot
    size=rep(1,length(chr))
    size[length>2000000]=2
    size[length>5000000]=3
    size[length>10000000]=4
    
    #Select the chromosome
    this <- chroms[chroms$c==c,]
    
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
         xlim = c(-1,1.5),
         ylim = ylim)
    
    #Insert Y and X axis
    axis(side=2,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=0,to=1,by=0.2),las=1)
    axis(side=1,cex.axis=0.6,tck=0.963,col.ticks='#808080',at=seq(from=-1,to=1.5,by=0.1))
    
    #Titles,date/time and axis labels
    mtext(text="log-ratio",side=1,line=1.1,cex=1)
    mtext(text="Allelic imbalance",side=2,line=1.5,cex=1)
    mtext(paste("Sample: ",name,"\n Chromosome ",c,", Region: ",Rstart,"-",Rend,sep=""),side=3)
    
    #------------------------------------------------------------
    #Top Right - Signal
    #------------------------------------------------------------
    
    #Go to screen 7
    screen(7)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(0.5,0,0))
    par(lend=1)
    par(xpd=T)
    
    #Select the correct chromosome and remove stuff lower than -1 
    mix <- mchr==as.character(this$chr) #& mval>(-1)
    
    #Create an index of colors relating to the positions and lengths on this chromosome
    col=rep('#000000',sum(ix))
    col[pos[ix] < this$mid] <- colors_p(sum(pos[ix] < this$mid))
    col[pos[ix] > this$mid] <- colors_q(sum(pos[ix] > this$mid))
    
    #Predefine ymin ymax and sequence between them
    ymin=floor(min(int[ix]))-0.5
    if(ymin > -1)
    {
        ymin = -1
    }
    ymax=ceiling(max(int[ix]))+0.5
    if(ymax < 1)
    {
        ymax = 1
    }
    seqminmax=seq(ymin,ymax,by=0.5)
    
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
         ylim = c(ymin,ymax)) #c(-1,1.5))
    
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
    mtext("log-ratio",side=2,line=0.3)
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
    #Bottom Right - Allele frequency
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
         labels=paste(seq(5,floor(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
    
    #Add a rectangle showing which region has been chosen
    rect(xleft=Rstart,xright=Rend,ybottom=0,ytop=1,col='#E8000040',lty=0)
    
    #Add X and Y label
    mtext("Allele frequency",side=2,line=0.3)
    mtext("Position (Mb)",side=1,line=1)
    
    #------------------------------------------------------------
    #Top Bottom - Detailed region signal view
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
    ymin=floor(min(int[wix]))-0.5
    if(ymin > -1)
    {
        ymin = -1
    }
    ymax=ceiling(max(int[wix]))+0.5
    if(ymax < 1)
    {
        ymax = 1
    }
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
    mtext("log-ratio",side=2,line=-0.8)
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
    #      labels=paste(seq(5,floor(this$length/1e6),5),'',sep=''),cex.axis=0.6,pos=0)
    axis(side=1,cex.axis=0.6)
    
    #Add X and Y label
    mtext("Allele frequency",side=2,line=-0.8)
    mtext("Position (Mb)",side=1,line=1)
    
    
    #Close all the opened split.screens and release the figure
    close.screen(all.screens=T)
    dev.off()
    setwd("..")
} , silent=F)
    
}

#Wrapper for getEstimates(). It will execute getEstimates in every sample folder.
#If SampleData.csv already exists it will create the file SampleData_YYYY-MM-DD_HH-MM-SS.csv instead.
TAPS_estimates <- function(path=getwd()) {
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(TAPS))
    # setwd("/media/safe/COAD_TCGA/NexusTxtFiles/tumorCELFiles")
    root <- path
    samples <- dir(path)[file.info(cbind(paste(path,cbind(dir(path)),sep='/')))$isdir]
    # samples <- dir(path)[file.info(path)$isdir]
    # samples <- samples[grep('^[A-Z]',samples)]
    datat <- foreach(sample=samples) %do% {
        setwd(paste(root,'/',sample,sep=''))
        if(file.exists('shortRegions.Rdata')) {
            load('shortRegions.Rdata')
            cn <- getEstimates(regs$logs,regs$scores)
            try(rm(regs),T)
            data.frame(Sample = sample,cn1=cn[1],cn2=cn[2],cn3=cn[3],loh=cn[4])
        } else {
            data.frame(Sample = sample,cn1=NA,cn2=NA,cn3=NA,loh=NA)
        }
    }
    sampledata <- do.call(rbind,datat)
    setwd(root)
    # if (file.exists('SampleData.csv')) {
    #     write.table(sampledata,paste('SampleData',format(Sys.time(), "_%F_%T.csv"),sep=''),sep='\t',row.names=F)
    # } else {
    #     write.table(sampledata,'SampleData.csv',sep='\t',row.names=F)
    # }
    if (file.exists('SampleData.csv')) {
            file.rename('SampleData.csv',paste('SampleData',format(Sys.time(), "_old_TAPS_estimates_%F_%T.csv"),sep=''))
    }
    write.table(sampledata,'SampleData.csv',sep='\t',row.names=F)
    print('OK')

}


getEstimates <- function(logR, imba, cellLines=F) {
    #load('shortRegions.Rdata')
    #logR=allRegions$regions$log2[allRegions$regions$lengthMB>5]
    #imba=allRegions$regions$imba[allRegions$regions$lengthMB>5]
    #logR=regs$logs
    #imba=regs$scores
    
    # Find the max point, that is likely an integer copy number.
    d=density(logR[abs(logR)<0.15],na.rm=T)
    max=d$x[order(d$y,decreasing=T)][1]
    
    # Find optimal logR's for integer copy numbers
    R=2^logR-2^max
    deltas=seq(0.05,0.33,0.01) # Allow a delta-ratio of 5% to 50%
    n=length(deltas)
    dev=rep(NA,n)
    for (i in 1:n) { 
        delta=deltas[i]
        mod=R %% delta
        mod[mod>(delta/2)]=delta-mod[mod>(delta/2)]
        dev[i]=mean(mod)
    }
    
    best=deltas[dev/deltas==min(dev/deltas)]
    
    # guess the copy number at "max"
    imbas=imba[abs(R)<best/4]
    imbas=imbas[!is.na(imbas)]
    d=density(imbas,na.rm=T)
    n=length(d$y)
    maxes=which(diff(d$y[1:(n-1)])>0 & diff(d$y[2:n])<0)
    ## remove crappy maxes:
    maxes=maxes[d$y[maxes] > 0.1*max(d$y[maxes])]
    
    cn=3
    if (d$x[maxes][1]<0.15) # even copy number
        cn=2
    if (length(maxes)>2)
        if (min(diff(maxes))/max(diff(maxes))>0.7)
            cn=4
    if (length(maxes)==1)
        if (d$x[maxes]>0.35)
            cn=1
    
    try(cn1 <- log2(2^max-(cn-1)*best),silent=T)
    try(cn1 <- median(logR[abs(logR-cn1)<0.1]),silent=T)
    cn2=log2(2^max-(cn-2)*best)
    try(cn2 <- median(logR[abs(logR-cn2)<0.1]),silent=T)
    cn3=log2(2^max-(cn-3)*best)
    try(cn3 <- median(logR[abs(logR-cn3)<0.1]),silent=T)

    
    R2=2^cn2
    imba2=imba[abs(2^logR-R2)<best/4]
    if(sum(is.na(imba2)) > (length(imba2)/2) | length(imba2) < 4) return(c(NA,NA,NA,NA))
    d=density(imba2,na.rm=T)
    n=length(d$y)
    maxes=which(diff(d$y[1:(n-1)])>0 & diff(d$y[2:n])<0)
    ## remove crappy maxes:
    maxes=maxes[d$y[maxes] > 0.05*max(d$y[maxes])]
    if (length(maxes)==1) {
        loh=0.75
    } else loh=d$x[maxes[length(maxes)]]
    return(round(c(cn1,cn2,cn3,loh),2))
}



#      _ _      _    
#  ___| (_) ___| | __
# / __| | |/ __| |/ /
#| (__| | | (__|   < 
# \___|_|_|\___|_|\_\
#                    
#click

drawImage <- function(file,dev=T,xlim,ylim) {
        jpeg <- readJPEG(file)
        if(dev==T) dev.new(height=11,width=15.5)
        par(mar=c(2,3,1,1),las=1,xaxs='i',yaxs='i',lend=1)
        plot(1, type="n", xlim=xlim, ylim=ylim,axes=F)
        rasterImage(jpeg,xlim[1], ylim[1], xlim[2], ylim[2], interpolate =T)
        # axis(1,at=seq(-3,3,0.25),lend=1)
        # axis(2,at=seq(ylim[1],ylim[2],0.05),lend=1)
}

clickWrap <- function(n,xlim,ylim) {
        y <- list()
        opar <- par
        for(i in 1:n) {
                print(i)
                col <- ifelse(i==1,'green','red')
                y[[i]] <- clicker(col) 
        }
        dev.new()
        par(mfrow=c(n,1),mar=c(2,2.5,.5,.5),xaxs='i',yaxs='i')
        lapply(y, FUN = function(y) {
               x <- 1:length(y)
               fit <- lm(y ~ x)
               xlim <- c(0.5,length(x)+.5)
               ylim <- c(min(y)-min(y)*.1,max(y)+max(y)*.1)
               plot(x,x*fit$coefficients[2]+fit$coefficients[1],type='l',ylim=ylim,xlim=xlim,las=1,lend=1)
               par(new=T)
               plot(x,y,xlab='',ylab='',axes=F,ylim=ylim,xlim=xlim,pch=20,cex=1.5)
               legend('topright',legend=paste('R^2: ',round(as.numeric(summary(fit)[8]),4),sep=''))
        })
        locator(1)
        try(dev.off(),T)
        par(opar)
}

clicker <- function(col='black') {
        print('Click on your values. Right click to close (This will redraw the image).')
        print(col)
        lastLogR <- locator(1,type='p',pch=20,col=col,cex=1.8)$x
        y <- 2^lastLogR
        cat(paste('#','logR','R','delta-R','\n',sep='\t'))
        cat(paste(1,round(lastLogR,2),round(2^lastLogR,2),0,'\n',sep="\t"))
        j <- 2
        while(T) {
                logR <- locator(1,type='p',pch=20,col=col,cex=1.8)$x
                if(is.null(logR)) break
                R <- 2^logR 
                y <- append(y,R)
                lastR <- 2^lastLogR
                diff <- R - lastR
                cat(paste(j,round(logR,3),round(R,3),round(diff,3),'\n',sep="\t"))
                lastLogR <- logR      
                j = j + 1
        }
        return(y)
}

copyNumberClicker <- function(maxy,miny) {
        coord <- list()
        for(i in 1:4) {
                if(i<4) {
                        print(paste("Left click on cn",i,". Right click if it doesn't exist",sep=''))
                } else {
                        print('Please select LOH')
                }
                tmp <- locator(n=1,type='p',pch=20,col='green',cex=1.8)
                if (is.null(tmp$x)) tmp$x <- tmp$y <- NA
                text(tmp, labels=c('cn1','cn2','cn3','LOH')[i],adj=1.5)
                text(tmp, labels=if(i < 4) round(tmp$x,2) else round((tmp$y - miny) / (maxy - miny),2) ,adj=c(1.3,+2.2))
                coord$x[i] <- tmp$x
                coord$y[i] <- tmp$y
        }
        return(coord)
}


returnSkippedValues <- function(sample,sampleDataOri,done,skipped){ 
    index <- sampleDataOri$Sample == sample
    return(list(sample,sampleDataOri[index,2],sampleDataOri[index,3],sampleDataOri[index,4],sampleDataOri[index,5],done,skipped))
}

TAPS_click <- function(path = getwd()) {

    suppressPackageStartupMessages(library(tcltk))
    suppressPackageStartupMessages(library(jpeg))
    suppressPackageStartupMessages(library(foreach))
    # require(xlsx)

    root <- path
    chr <- paste('Chromosome',c(1:22,'X'))
    xlim <- c(-2.387, 2.16)
    ylim <- c(0, 1)
    maxy <- 0.9336337
    miny <- 0.5170546

    # sampleData <- read.csv('SampleData.csv',sep='\t',header=T,colClasses=c('character',rep('numeric',4)),stringsAsFactors=F)
    sampleData <- read.csv(paste(path,'SampleData.csv',sep='/'),sep='\t',header=T,stringsAsFactors=F)
    if(any(colnames(sampleData) %in% c('MAPD','MHOF'))) sampleData <- sampleData[,1:5]
    sampleDataOri <- sampleData

    if(any(T == sampleData$done)) {
        print(head(sampleData,25))
        ans <- menu(c('Continue','Restart'),title='\n"Done"-column detected.\nWould you like to start from the beginning or continue from where you left off?')
        if(ans==1) {
            sampleData <- sampleData[sampleData$done == F, ]
            newSampleData <- sampleData
            restart <- F
        } else {
            newSampleData <- sampleData
            newSampleData$done <- F
            newSampleData$skipped <- F
            restart <- T
        }
    } else {
        newSampleData <- sampleData
        newSampleData$done <- F
        newSampleData$skipped <- F
        restart <- F
    }
        

    #index <- sampleData$Sample %in% examinedSampleData$Sample[examinedSampleData$done == T]
    #sampleData <- sampleData[!index,]
    #sampleData <- examinedSampleData[apply(examinedSampleData[,2:4],1,FUN=function(x) {sum(is.na(x))}) == 2,]

    setwd(root)
    close <- F

    write.table(x=rbind(colnames(newSampleData)),file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,col.names=F)
    if(nrow(sampleDataOri[sampleDataOri$done==T,]) > 0) {
        write.table(x=sampleDataOri[sampleDataOri$done==T,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,col.names=F,append=T)
    }

    copyNumbers <- foreach(sample=sampleData[,1] ,cn1=sampleData[,2],cn2=sampleData[,3],cn3=sampleData[,4],loh=sampleData[,5]) %do% {
        if(close == T) {

            # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,F,T)
            newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,T)

            write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
            return(newSampleData[newSampleData$Sample == sample,])
        }
        print(sample)
        setwd(paste(root,sample,sep='/'))
        file=dir()[grep('chr1.jpg',dir())]

        # if(sum(is.na(c(cn1,cn2,cn3))) > 1 ) {
        if(is.na(loh)) {
            drawImage(file,T,xlim=xlim,ylim=ylim)
        }

        skip <- F
        i <- 1
        lastAnswer <- 0
        coord <- list()
        selectCn <- F

        while(i > 0) {
            # if(sum(is.na(c(cn1,cn2,cn3))) > 1) { #If no preexisting value exist 
            # browser()
            if(is.na(loh)) { #If no preexisting value exist 
                if(lastAnswer == 10) {
                    drawImage(file,xlim=xlim,ylim=ylim)
                    lastAnswer <- 0
                }

                # Skip the menu and go directly to the copy number choice
                if(selectCn == T) {
                    answer <- 1
                    selectCn = F
                } else {
                    answer <- menu(c('Choose copynumbers','Sample overview',chr,'Check separation','Skip sample',"Close program",'Browser()'),
                                   graphics=TRUE,title=paste(sample,' (', which(sampleData$Sample==sample) ,'/' ,dim(sampleData)[1],')',sep=''))
                }

                tclServiceMode(TRUE)

                if(answer == 1) { #Copy number menu
                    if(lastAnswer == 2) drawImage(file,F,xlim=xlim,ylim=ylim)
                    coord <- copyNumberClicker(maxy,miny)
                    dev.off()
                    i <- 10
                } else if(answer == 2) { #Overview menu
                    overview=dir()[grep('overview.jpg',dir())]
                    drawImage(overview,F,xlim=xlim,ylim=ylim)
                    if(lastAnswer == 0) lastAnswer <- 2
                } else if(answer >= 3 & answer <=24) { #Select chromosome #1-22
                    file=dir()[grep(paste('chr',answer-2,'.jpg',sep=''),dir())]
                    drawImage(file,F,xlim=xlim,ylim=ylim)
                } else if(answer == 25) { #Select the X-chromosome
                    file=dir()[grep(paste('chr','X','.jpg',sep=''),dir())]
                    drawImage(file,F,xlim=xlim,ylim=ylim)
                } else if(answer == 26) { #Check separation
                    ans <- menu(c('1       ','2    '),title='1 or 2?',graphics=T)
                    tclServiceMode(TRUE)
                    clickWrap(ans,xlim=xlim,ylim=ylim)              
                    lastAnswer <- 10
                    try(dev.off(),T)
                } else if(answer == 27 | answer == 0) { #Skip sample
                    try(dev.off(),T)
                    # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,T,T)
                    newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,T)
                    write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
                    return(newSampleData[newSampleData$Sample == sample,])
                } else if(answer == 28) { #Close everything!
                    close <- T
                    skip <- T
                    i <- 0
                    try(dev.off(),T)
                    # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,F,T)
                    newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,T)
                    write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
                    return(newSampleData[newSampleData$Sample == sample,])
                } else if(answer == 29) { #Get a browser()
                    #browser()
                }

            } else { #If pre-existing values for the sample exists, grab the values.
                i <- 10
                lohCoord <- loh * (maxy - miny) + miny
                coord <- list(x=c(cn1,cn2,cn3,0),y=c(0,0,0,lohCoord) )
                cn1 <- cn2 <- cn3 <- loh <- NA
            }

            if(i == 10) { # i will be equal to 10 when the user have choosen the copy numbers and will be promted if the values are correct.
                drawImage(file,xlim=xlim,ylim=ylim)
                LOH <- (coord$y[4] - miny) / (maxy - miny)
                text(x=coord$x[1:3]-0.1,y=miny+0.016,labels=round(c(coord$x[1:3]),2))
                text(x=coord$x[1:3]-.1,y=miny+0.036, c('cn1','cn2','cn3'))

                #CN
                segments(coord$x[1:3],miny,coord$x[1:3],maxy,lwd=1.5,lend=1)
                segments(-2,coord$y[4],2,coord$y[4],lwd=1.5,lend=1)

                #LOH
                text(x=-1.82,y=coord$y[4]+.012,labels=round(c(LOH),2),cex=1,1)
                text(x=-1.82,y=coord$y[4]+.032, c('LOH'),cex=1,1)

                answer2 <- menu(c('Next sample','Redo sample','Redo copy number selection','Skip sample','Close program'),
                                title=paste(sample,' (', which(sampleData$Sample==sample) ,'/' ,dim(sampleData)[1],')',sep=''),graphics=T)

                if (answer2 == 1) {            # Next sample
                    i <- 0
                    size <- dev.size(units='px')
                    dev.copy(png,filename=paste(sample,'_selected_cn_and_LOH.png',sep=''),width=size[1],height=size[2],units='px')
                    dev.off()

                } else if (answer2 == 2) {     #Redo sample
                    i <- 1
                    print(sample)
                } else if (answer2 == 3) {     #Redo copy number selection
                    i <- 1
                    selectCn <- T
                    print(sample)
                } else if (answer2 == 4) {     #Skip sample
                    try(dev.off(),T)
                    # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,T,T)
                    newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,T)
                    write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
                    return(newSampleData[newSampleData$Sample == sample,])
                } else if (answer2 == 5) {     #Close program
                    close <- T
                    skip <- T
                    i <- 0
                    try(dev.off(),T)
                    # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,F,T)
                    newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,T)
                    write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
                    return(newSampleData[newSampleData$Sample == sample,])
                }
                lastAnswer <- 10
                try(dev.off(),silent=T)
            }
        }
        # Write data to the data.frame
        if(skip==F) {
            LOH <- (coord$y[4] - miny) / (maxy - miny)
            newSampleData[newSampleData$Sample == sample,] = list(sample,coord$x[1],coord$x[2],coord$x[3],LOH,T,skip)
            write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
            return(newSampleData[newSampleData$Sample == sample,])
        } else {
            # newSampleData[newSampleData$Sample == sample,] = list(sample,NA,NA,NA,NA,T,skip)
            newSampleData[newSampleData$Sample == sample,] = returnSkippedValues(sample,sampleDataOri,F,skip)
            write.table(x=newSampleData[newSampleData$Sample == sample,],file=paste(root,'backup.csv',sep='/'),sep='\t',row.names=F,append=T,col.names=F)
            return(newSampleData[newSampleData$Sample == sample,])
        }
        tclServiceMode(TRUE)
        try(dev.off(),silent=T)
        1
    }
    setwd(root)
    if (file.exists('SampleData.csv')) {
        file.rename('SampleData.csv',paste('SampleData',format(Sys.time(), "_old_TAPS_click_%F_%T.csv"),sep=''))
    }
    append <- nrow(sampleDataOri[sampleDataOri$done==T,])>0

    if(restart == F) { 
        write.table(sampleDataOri[sampleDataOri$done==T,],'SampleData.csv',sep='\t',row.names=F)
    } 

    if (restart == T ) {
        write.table(do.call(rbind,copyNumbers), 'SampleData.csv',sep='\t',row.names=F,append=F,col.names= T)
    } else {
        write.table(do.call(rbind,copyNumbers), 'SampleData.csv',sep='\t',row.names=F,append=append,col.names= !append)
    }

}

chromFix <- function(data) {
    if (any(data$Chromosome=='1'))
        data$Chromosome=paste('chr',data$Chromosome,sep='')
    return(data)
}

