patchwork.copynumbers = function(CNfile,cn2,delta,het,hom,maxCn=8,ceiling=1,forcedelta=F,male.sample=F,male2femref=F)
	{
	
	#packagepath = system.file(package="patchwork")
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))
	#load(paste(packagepath,"/data/ideogram.RData",sep=""))
	data(ideogram,package="patchworkData")
	load(CNfile)

    #save parameters as strings
    parameters=paste("Parameters given: cn2:",cn2," delta:",delta," het:",het," hom:",hom)

    name = strsplit(tolower(CNfile),"_copynumbers.rdata")

    name = strsplit(name[[1]],"/")

    name = name[[1]][length(name[[1]])]

	voidchrom <- c('chrX','chrY') # may add non-integer chroms here....

    tix=NULL     #temporary index
    int=NULL     ## contains coverage estimate of each (total) copy number
    ai=NULL         ## contains Allelic Imbalance Ratio estimate of each copy number variant.
    regions <- segs
    regions <- regions[(is.autosome(regions$chr)&regions$np>50)&(!is.nan(regions$ai)),] ## will use these regions
    regions <- regions[!is.na(regions$median),]

    ## likely cn2 regions sit within delta/3 from cn2.
    expectedAt <- cn2
    tix$cn2 <- abs(regions$median - expectedAt) < (delta/3)    ## index of likely cn2 regions
    temp <- regions[tix$cn2,]                                ## cn2 regions
    med <- weightedMedian(temp$median,temp$np)            ## improved value of Log-R at cn2 (returns NULL if theres nothing there)
    int$cn2 <- ifelse(!is.null(med),med,expectedAt)            ## saved to int.

    ## likely cn1 regions sit at about cn2 - delta:
    d <- delta                                                 ## the (Log-R) distance to cn1
    expectedAt <- int$cn2-d                                     ## cn1 is expected here
    tix$cn1 <- abs(regions$median - expectedAt) < (d/3)        ## index of likely cn1 regions
    temp <- regions[tix$cn1,]
    med <- weightedMedian(temp$median,temp$np)
    int$cn1 <- ifelse(!is.null(med),med,expectedAt)
    if (forcedelta) int$cn1 <- int$cn2-delta

    ## likely cn0 regions sit below cn1 - delta:
    d <- int$cn2-int$cn1
    expectedAt <- int$cn1-d                                     
    tix$cn0 <- regions$median < expectedAt
    temp <- regions[tix$cn0,]
    med <- weightedMedian(temp$median,temp$np)
    int$cn0 <- ifelse(!is.null(med),med,expectedAt)
    if (forcedelta) int$cn0 <- int$cn1-delta


    ## likely cn3 regions sit at about cn2+delta
    d <- delta
    expectedAt <- int$cn2+d
    tix$cn3 <- abs(regions$median - expectedAt) < (d/3)
    temp <- regions[tix$cn3,]
    med <- weightedMedian(temp$median,temp$np)
    int$cn3 <- ifelse(!is.null(med),med,expectedAt)
    if (forcedelta) int$cn3 <- int$cn2+delta


    ## cn4 follows at ...
    d <- delta
    expectedAt <- int$cn3+d
    tix$cn4 <- abs(regions$median - expectedAt) < (d/4)
    temp <- regions[tix$cn4,]
    med <- weightedMedian(temp$median,temp$np)
    int$cn4 <- ifelse(!is.null(med),med,expectedAt)
    if (forcedelta) int$cn4 <- int$cn3+delta

    ## generalized for higher cns
    for (cn in 5:maxCn)
    	{
    	thisCn <- paste('cn',cn,sep='')
    	prevCn <- paste('cn',cn-1,sep='')
    	pprevCn <- paste('cn',cn-2,sep='')
    	d <- delta #dmin*(int[prevCn][[1]]-int[pprevCn][[1]])
    	expectedAt <- int[prevCn][[1]]+d
    	tix[[thisCn]] <- abs(regions$median - expectedAt) < (d/5)
    	temp <- regions[tix[thisCn][[1]],]
    	med <- weightedMedian(temp$median,temp$np)
    	int[thisCn] <- ifelse(!is.null(med),med,expectedAt)
        if (forcedelta) int[[thisCn]] <- int[[prevCn]]+delta
    	}

	#  smooth the cn relationship ?

    ## at cn2, find the variant clusters (normal and CNNLOH)
    ix <- (abs(regions$median - int$cn2) < 0.2*(int$cn3-int$cn2) ) # taking only closely-matching segments
    data <- regions[ix,]
    data <- data[!is.na(data$ai),] # ...with a calculated allelic imbalance.

    expectedAt <- het
    ix <- 2*abs(data$ai-het) < abs(data$ai-hom)
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    ai$cn2m1 <- ifelse (!is.null(med),med,expectedAt)

    expectedAt <- hom
    ix <- 2*abs(data$ai-hom) < abs(data$ai-het)
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    ai$cn2m0 <- ifelse (!is.null(med),med,expectedAt)

    ## for cn1 (and 0)
    ix <- (abs(regions$median - int$cn1) < 0.2*(int$cn2-int$cn1) )
    data <- regions[ix,]
    data <- data[!is.na(data$ai),]
    expectedAt <- (ai$cn2m0+ai$cn2m1)*3/5     ## Decent estimate.
    med <- weightedMedian(data$ai,data$snvs) ## Average allelic imbalance weighted on snp count
    ai$cn1m0 <- ifelse (!is.null(med),med,expectedAt)    ## Will be NA if there was no CNNLOH
    ai$cn0m0 <- NA #unimportant

    ## for cn3:
    ix <- (abs(regions$median - int$cn3) < 0.2*(int$cn4-int$cn3) )
    data <- regions[ix,]
    data <- data[!is.na(data$ai),]
  
    range <- ai$cn2m0-ai$cn2m1 #the distance between 2normal and CNNLOH
    # get the 3(1) regions: 
    expectedAt <- ai$cn2m1+range/3 # this is an approx if cn3m1 is absent.

    ix <- (data$ai<ai$cn2m0) & (data$ai>ai$cn2m1) # take regions with less AI than cn2m0 but more than cn2m1
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) # average allelic imbalance weighted on snp count
    ai$cn3m1 <- ifelse (!is.null(med),med,expectedAt)

    # now for cn3m0
    expectedAt <- ai$cn3m1 + range # approx of cn3m0
    ix <- (abs(data$ai-expectedAt) < (expectedAt-ai$cn3m1)/4)  # take those much closer to exp than to cn3m1
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    try (if (med<ai$cn2m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
    ai$cn3m0 <- ifelse (!is.null(med),med,expectedAt)

    if (is.na(ai$cn1m0)) ai$cn1m0 <- mean(ai$cn3m1,ai$cn2m0) ## If deletions were missing, place an estimate from cn3 
    
    ## now for cn4
    ix <- abs(regions$median - int$cn4) < 0.2*(int$cn4-int$cn3) 
    data <- regions[ix,]
    data <- data[!is.na(data$ai),]

    # get the 4(2) regions: they are at AI about cn2m1
    expectedAt <- ai$cn2m1 # this is a good approx 

    ix <- data$ai<ai$cn3m1 # let all below cn3m1 in
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    ai$cn4m2 <- ifelse (!is.null(med),med,expectedAt)
    if (ai$cn4m2>ai$cn2m1) ai$cn4m2 <- ai$cn2m1 # sanity check.

    data <- data[!ix,] # coutinue with the remaining
    # now 4(1) has less ai than 3(0) -> sits at about cn3m1+(cn3m1-cn2m1).
    expectedAt <- ai$cn3m1+(ai$cn3m1-ai$cn2m1)
    ix <- abs(data$ai-expectedAt)<abs(data$ai-ai$cn3m0) # take those closer to exp than to 3(0)
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    ai$cn4m1 <- ifelse (!is.null(med),med,expectedAt)

    # now 4(0) is at about 3(0) + [0.9 - 3(0)] * [3(0)-2(0)] / [0.9-2(0)] 
    expectedAt <- ai$cn3m0 + (ceiling-ai$cn3m0)*(ai$cn3m0-ai$cn2m0)/(ceiling-ai$cn2m0)
    ix <- abs(data$ai-expectedAt) < 0.2*(expectedAt-ai$cn4m2) # we take those very close to exp
    med <- weightedMedian(data$ai[ix],data$snvs[ix])
     try (if (med<ai$cn3m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
    ai$cn4m0 <- ifelse (!is.null(med),med,expectedAt)

    ## generalization for higher copy numbers. it gets complicated.
    for (cn in 5:maxCn) {
    thisCn <- paste('cn',cn,sep='')
    prevCn <- paste('cn',cn-1,sep='')
    pprevCn <- paste('cn',cn-2,sep='')
    
    ix <- (abs(regions$median - int[thisCn][[1]]) < 0.2*(int[thisCn][[1]]-int[prevCn][[1]]) )
    data <- regions[ix,]
    data <- data[!is.na(data$ai),]
    data <- data[data$np>50,] # long regions for safety

    ## try to find variants, starting with LOH
    # LOH such as 5(0)
    m <- 0
    thisVariant=paste(thisCn,'m',0,sep='')
    c4m0 <- ai[paste(prevCn,'m',m,sep='')][[1]] # relative naming for clarity
    c3m0 <- ai[paste(pprevCn,'m',m,sep='')][[1]]
    
    expectedAt <- ceiling-((ceiling-c4m0)*(ceiling-max(c4m0,c3m0))/(ceiling-min(c3m0,c4m0)))
    #ai[thisVariant] <- expectedAt
    # we then take those close to exp (within delta[3(0),4(0)])
    ix <- abs(data$ai-expectedAt) < abs(c4m0 - c3m0) 
    med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
    try (if (med<ai$cn4m0) med <- NULL, silent=T)    ## in case of heterogeneities or strange signals, we disallow this effect.
    ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
    
    ## then from balanced to less balanced
    minorVariants=trunc(cn/2):1 
    first <- T
    for (m in minorVariants) {
        thisVariant=paste(thisCn,'m',m,sep='')
        if (m==cn/2) {
        # We have balanced variant, so rather easy.
        expectedAt <- ai[paste(pprevCn,'m',m-1,sep='')][[1]] # this is a good approx, balanced at cn-2
        ix <- data$ai<ai[paste(prevCn,'m',m-1,sep='')][[1]] # let all below (cn-1, mcn-1) in
        med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
	if (ai[thisVariant] > ai[paste(pprevCn,'m',m-1,sep='')][[1]]) ai[thisVariant] <- ai[paste(pprevCn,'m',m-1,sep='')][[1]] # sanity check.
        first <- F
        } else if (first) { 
        # its not balanced but its the most balanced of the unbalanced. something like 5(2)
        expectedAt <- 0.5*( ai[paste(prevCn,'m',m,sep='')][[1]] + ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) # that means between 4(2) and 3(1)
        ix <- abs(data$ai-expectedAt) < ( ai[paste(prevCn,'m',m,sep='')][[1]] - ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) /3 # let all "between" 4(2) and 3(1) in
        med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
        first <- F
        } else {
        # not the most balanced unbalanced variant, for example 5(1):
        expectedAt <- ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]] + (minorVariants[1]-m) * (ai[paste(thisCn,'m',0,sep='')][[1]] - ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]]) / trunc(cn/2) # 5(2) + (which)* 5(0)-5(2) /(n)
        ix <- abs(data$ai-expectedAt) < ( ai[paste(prevCn,'m',m,sep='')][[1]] - ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) /3 # let all "between" 4(1) and 3(0) in
        med <- weightedMedian(data$ai[ix],data$snvs[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
        } 
    } # done with minor variants
    } # done with copy numbers

    sampleInfo <- list('int'=int,'ai'=ai)

	##  Assign total and minor copy numbers to all segments.
    regions <- segs[!is.na(segs$median),]    ## This time, work on all segments available.

    Cn <- NULL            ## Total copy number
    mCn <- NULL            ## Minor allele copy number
    fullCN <- NULL        ## Variant label. ('cnXmY')

    intDist <- NULL        ## distance to certain Log-R
    imbaDist <- NULL    ## distance to certain allelic imbalance

	

    ## give each segment the best matching CN and mCN
    for (i in 1:nrow(regions)) 
    	{
   	 	# set total copy number
    	distance <- Inf
    	for (cn in 0:maxCn) 
    		{
        	t_int <- int[paste('cn',cn,sep='')][[1]]    ## get Log-R of particular cn from 'int'
        	t_dis <- abs(regions$median[i]-t_int)            ## distance to that particular cn
        	if (t_dis < distance) 
        		{                        ## nearest so far, save.
        		distance <- t_dis -> intDist[i]
        		Cn[i] <- cn
        		}
    		}
    	}
    # Y makes CN0 sit very low, this is a fix on non-Y Cn0
    #if ((Cn[i] == 1)&(intDist[i] > (int$cn2-int$cn1))) Cn[i] <- 0 ## currently not needed.
	for (i in 1:nrow(regions)) 
		{
   		# set minor CN
    	distance <- Inf
    	if (Cn[i]<=1)
    		{
        	mCn[i] <- 0
    		} 
    	else if (!is.na(regions$ai[i])) 
    		if (!is.nan(regions$ai[i])) 
    			for (m in 0:trunc(Cn[i]/2)) 
    				{
        			t_ai <- ai[paste('cn',Cn[i],'m',m,sep='')][[1]]
        			t_dis <- abs(regions$ai[i]-t_ai)
        			if (t_dis < distance) 
        				{
        				distance <- t_dis -> imbaDist[i]
        				mCn[i] <- m
        				}
    				}
    		else 
    			mCn[i] <- NA
    			fullCN[i] <- paste('cn',Cn[i],'m',mCn[i],sep='') # Full description
  		}
 

    regions$Cn <- Cn
    regions$mCn <- mCn
    regions$fullCN <- fullCN

    temp <- regions[is.autosome(regions$chr),]
    temp <- temp[!is.na(temp$Cn),]
	
 
    meanCn <- weightedMean(temp$Cn,temp$np)
	ix1 <- temp$Cn==trunc(meanCn)
	ix2 <- temp$Cn==ceiling(meanCn)
	xdelta <- weightedMean(temp$median[ix2],temp$np[ix2]) - weightedMean(temp$median[ix1],temp$np[ix1])
    xdelta <- xdelta / weightedMean(temp$median,temp$np)

	expected_delta <- 1/meanCn
	
	tumor_percentDNA <- xdelta / expected_delta
    tumor_percent <- tumor_percentDNA/(meanCn/2) / ( tumor_percentDNA/(meanCn/2) + (1-tumor_percentDNA) )

    ## Fix sex chromosomes in case of male sample.
    if (male.sample) {
        ix <- which(!is.autosome(regions$chr))
        regions$mCn[ix] <- 0
        ## Treats differently depending on sex of the reference
        if (!male2femref) {
            # X and Y (originally cn1) were compared to cn1 (matched normal or male references).
            # Unchanged copy number (=1) is expected at relative coverage 1. 
            # Gains/loss are expected at xdelta (the delta of the autosomes)*2.
            rawcopychange <- (regions$median[ix]-1) / (xdelta*2)
            regions$Cn[ix] <- round(1+rawcopychange)
            regions$fullCN[ix] <- paste('cn',round(1+rawcopychange),'m0')
        } else 
        {
            # X and Y (originally cn1) were compared to cn2 (female references).
            # Unchanged copy number (=1) is expected at relative coverage 1/2. 
            # Gains/loss are expected at xdelta from there.
            rawcopychange <- (regions$median[ix]-1/2) / (xdelta)
            regions$Cn[ix] <- round(1+rawcopychange)
            regions$fullCN[ix] <- paste('cn',round(1+rawcopychange),'m0')
        }
    }
    ## And a small fix on female samples
    if (!male.sample) {
        ix <- which(regions$chr=='chrY')
        regions$mCn[ix] <- regions$Cn[ix] <- 0
    }


    regions$meanCn <- meanCn
    regions$tumor_percent <- tumor_percent
    write.csv(regions,file=paste(name,'_Copynumbers.csv',sep=""))

    #make t just for intuitive work between patchwork and taps
    t=list(int=int,ai=ai)

    save(regions,kbsegs,alf,t,file="Development.Rdata")


	karyotype_check(regions$chr,regions$start,regions$end,regions$median,regions$ai,
					regions$Cn,regions$mCn,t,name=name,
					xlim=c(0,2.4),ylim=c(0.1,1))
					
	# karyotype_chromsCN(regions$chr,regions$start,regions$end,regions$median,regions$ai,
	# 					regions$Cn,regions$mCn,kbsegs$chr,
	# 					kbsegs$pos,kbsegs$ratio,alf$achr,alf$apos,
	# 					(1-alf$amin/alf$amax),name=name,xlim=c(0,2.4),
	# 					ylim=c(0.1,1))
    karyotype_chromsCN(regions$chr,regions$start,regions$end,regions$median,regions$ai,regions$Cn,regions$mCn,
        kbsegs$chr,kbsegs$pos,kbsegs$ratio,
        alf$achr,alf$apos,(1-alf$amin/alf$amax),t,
        name=name,parameters)
	}