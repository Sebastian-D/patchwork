patchwork.CG.plot <- function(path=NULL,name='CG_sample',manual_file_input = FALSE,masterVarBeta=NULL,somaticCnvSegments=NULL,depthOfCoverage=NULL)
	{
	
	#------------------------------------------------#
	# Read input files
	#------------------------------------------------#
	
	#Choose one of the three modes of input.
	if(manual_file_input==TRUE)
		{
		#Prompt for input of path to your files.
		cat("ex: ~/user/documents/projekt/masterVarBeta-123-456-ASM-T1-N1.tsv.bz2")
		cat("\n","Enter complete path AND filename to your masterVarBeta file:","\n")
		mv = scan(n=1,what=character())
		cat("\n","ex: ~/user/documents/projekt/somaticCnvSegmentsNondiploidBeta-123-456-ASM-T1-N1.tsv")
		cat("\n","Enter complete path AND filename to your somaticCnvSegmentsNondiploidBeta file:","\n")
		sm = scan(n=1,what=character())
		cat("\n","ex: ~/user/documents/projekt/depthOfCoverage_100000-123-456-ASM-T1.tsv")
		cat("\n","Enter complete path AND filename to your depthOfCoverage file:","\n")
		dc = scan(n=1,what=character())
		

		#Check if the files exists
		if(file.exists(mv)==F)
			{
			stop("The required file,masterVarBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(file.exists(sm)==F)
			{
			stop("The required file,somaticCnvSegmentsNondiploidBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(file.exists(dc)==F)
			{
			stop("The required file,depthOfCoverage, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}

		skip_mv = scan(file=mv,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">locus",skip_mv)
		while(length(j) != 1)
			{
			skip_mv = scan(file=mv,n=1,skip=i,what=character(),quiet=T)
			j = grep(">locus",skip_mv)
			i = i + 1
			}
		skip_mv = i + 1
		
		skip_sm = scan(file=sm,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_sm)
		while(length(j) != 1)
			{
			skip_sm = scan(file=sm,n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_sm)
			i = i + 1
			}
		skip_sm = i + 1
		
		skip_dc = scan(file=dc,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_dc)
		while(length(j) != 1)
			{
			skip_sm = scan(file=dc,n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_dc)
			i = i + 1
			}
		skip_dc = i + 1
		
		cat("\n","Reading files","\n")
		
		#Read masterVarBeta
		mastervar = read.table(mv, sep="\t",skip=skip_mv,
						colClasses=c("NULL","NULL","character","integer","integer",
						"NULL","character",rep("NULL",11),"character",rep("NULL",4),"numeric","numeric",
						rep("NULL",17),"integer","integer",rep("NULL",5)))
						
		#Read somaticCnvSegmentsNondiploidBeta
		segm = read.table(sm,sep="\t",skip=skip_sm,
						colClasses=c("character","integer","integer","numeric","numeric",
						rep("NULL",7)))
		
		#Read depthOfCoverage
		depcov = read.table(dc,sep="\t",skip=skip_dc,
							colClasses=c("character","integer","integer",rep("NULL",3),"numeric"))
				
		}
	else if(length(masterVarBeta) == 1 & length(somaticCnvSegments) == 1 & length(depthOfCoverage) == 1)
		{

		#Check if the files exists
		if(file.exists(masterVarBeta)==F)
			{
			stop("The required file,masterVarBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(file.exists(masterVarBeta)==F)
			{
			stop("The required file,somaticCnvSegmentsNondiploidBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(file.exists(masterVarBeta)==F)
			{
			stop("The required file,depthOfCoverage, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}

		mv_ = strsplit(masterVarBeta,"/")
		sm_ = strsplit(somaticCnvSegments,"/")
		dc_ = strsplit(depthOfCoverage,"/")
		
		mv = grep("masterVar",mv_[[1]],value=T)
		sm = grep("somaticCnv",sm_[[1]],value=T)
		dc = grep("depthOfCoverage",dc_[[1]],value=T)
		
		skip_mv = scan(file=masterVarBeta,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">locus",skip_mv)
		while(length(j) != 1)
			{
			skip_mv = scan(file=masterVarBeta,n=1,skip=i,what=character(),quiet=T)
			j = grep(">locus",skip_mv)
			i = i + 1
			}
		skip_mv = i + 1
		
		skip_sm = scan(file=somaticCnvSegments,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_sm)
		while(length(j) != 1)
			{
			skip_sm = scan(file=somaticCnvSegments,n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_sm)
			i = i + 1
			}
		skip_sm = i + 1
		
		skip_dc = scan(file=depthOfCoverage,n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_dc)
		while(length(j) != 1)
			{
			skip_dc = scan(file=depthOfCoverage,n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_dc)
			i = i + 1
			}
		skip_dc = i + 1
		
		cat("Reading files: \n",mv,"\n",sm,"\n",dc,"\n")
		#Read mastervarbeta file in ASM folder
		mastervar = read.table(masterVarBeta,
						sep="\t",
						colClasses=c("NULL","NULL","character","integer","integer",
						"NULL","character",rep("NULL",11),"character",rep("NULL",4),"numeric","numeric",
						rep("NULL",17),"integer","integer",rep("NULL",5)),
						skip=skip_mv)
						
		#Read somaticCnvSegmentsNondiploidBeta in ASM/CNV folder
		segm = read.table(somaticCnvSegments,
					sep="\t",skip=skip_sm,
					colClasses=c("character","integer","integer","numeric","numeric",
					rep("NULL",7)))
					
		depcov = read.table(depthOfCoverage,sep="\t",skip=skip_dc,
							colClasses=c("character","integer","integer",rep("NULL",3),"numeric"))
		}
	else
		{
		#Check that path is given
		if(length(path)==0)
			{
			stop("You must input the required files, either by the path parameter, manual_file_input parameter or individual pointer parameters to the files in question. See ?patchwork.CG.plot")
			}
		mastervar_ = list.files(path=path,pattern="masterVar")
		path_segm = paste(path,"/CNV",sep="")
		segm_ = list.files(path=path_segm,pattern="somaticCnvSegmentsNondiploid")
		depcov_ = list.files(path=path_segm,pattern="depthOfCoverage")
		
		#If the files were not found.
		if(length(mastervar_)==0)
			{
			stop("The required file,masterVarBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(length(segm_)==0)
			{
			stop("The required file,somaticCnvSegmentsNondiploidBeta, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}
		if(length(depcov_)==0)
			{
			stop("The required file,depthOfCoverage, could not be found. Did you input the correct path? See ?patchwork.CG.plot")
			}

		skip_mv = scan(file=paste(path,"/",mastervar_,sep=""),n=1,what=character(),quiet=T)
		i = 1
		j = grep(">locus",skip_mv)
		while(length(j) != 1)
			{
			skip_mv = scan(file=paste(path,"/",mastervar_,sep=""),n=1,skip=i,what=character(),quiet=T)
			j = grep(">locus",skip_mv)
			i = i + 1
			}
		skip_mv = i + 1
		
		skip_sm = scan(file=paste(path_segm,"/",segm_,sep=""),n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_sm)
		while(length(j) != 1)
			{
			skip_sm = scan(file=paste(path_segm,"/",segm_,sep=""),n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_sm)
			i = i + 1
			}
		skip_sm = i + 1
		
		skip_dc = scan(file=paste(path_segm,"/",depcov_,sep=""),n=1,what=character(),quiet=T)
		i = 1
		j = grep(">chr",skip_dc)
		while(length(j) != 1)
			{
			skip_dc = scan(file=paste(path_segm,"/",depcov_,sep=""),n=1,skip=i,what=character(),quiet=T)
			j = grep(">chr",skip_dc)
			i = i + 1
			}
		skip_dc = i + 1
	
	
		cat("Reading files from ASM folder \n")
	
		#Read mastervarbeta file in ASM folder
		mastervar = read.table(paste(path,"/",mastervar_,sep=""),
							sep="\t",
							colClasses=c("NULL","NULL","character","integer","integer",
							"NULL","character",rep("NULL",11),"character",rep("NULL",4),"numeric","numeric",
							rep("NULL",17),"integer","integer",rep("NULL",5)),
							skip=skip_mv)
	
		#Read somaticCnvSegmentsNondiploidBeta in ASM/CNV folder
		segm = read.table(paste(path_segm,"/",segm_,sep=""),
					sep="\t",skip=skip_sm,
					colClasses=c("character","integer","integer","numeric","numeric",
					rep("NULL",7)))
					
		depcov = read.table(paste(path_segm,"/",depcov_,sep=""),sep="\t",skip=skip_dc,
							colClasses=c("character","integer","integer",rep("NULL",3),"numeric"))
	
		}
		
		colnames(segm)=c("chr","start","end","avgnormcov","relcov")
							
		colnames(mastervar)=c("chr","begin","end","vartype","dbSNP","ref_count","tot_count",
									"ref_countN","tot_countN")
		
		colnames(depcov)=c("chr","begin","end","avgnormcov")
		
		cat("File input complete \n")
	
	#------------------------------------------------#
	# Filter out unecessary/unwanted values
	#------------------------------------------------#
	
	#Remove all rows that have NA values in column 9.(tot_countN).
	mastervar = mastervar[!is.na(mastervar[,9]),]
	
	#Only keep when referenceallelereadcount/totalreadcount is between 0.2 and 0.8
	ix = ((mastervar$ref_countN/mastervar$tot_countN) > 0.2 &
		(mastervar$ref_countN/mastervar$tot_countN) < 0.8)

	mastervar = mastervar[ix,]

	#Remove NA values generated by for example 0/0
	mastervar = mastervar[!is.na(mastervar[,9]),]

	#only keep markers wiht dbSNP id present.

	ix = grep("dbsnp",mastervar$dbSNP)

	mastervar = mastervar[ix,]
	
	#------------------------------------------------#
	# Add and calculate some columns for mastervar
	#------------------------------------------------#
	
	cat("Performing calculations \n")
	
	depcov$ratio <- NA
	
	depcov$ratio = depcov$avgnormcov / mean(depcov$avgnormcov)
	
	mastervar$ratio <- mastervar$min <- mastervar$max <- mastervar$mut_count <- NA

	#Calculate and fill min/max columns
	mastervar$mut_count = mastervar$tot_count - mastervar$ref_count
	mastervar$max = apply(mastervar[,c(6,10)],1,max)
	mastervar$min = mastervar$tot_count - mastervar$max
	mastervar$ratio = mastervar$tot_count / mean(mastervar$tot_count)
	
	# segm$ai <- segm$snvs <- segm$ratio <- NA

	# #Generate the allelic imbalance calculation, among other things.
	# for (i in 1:nrow(segm))
	# 	{
	# 	#Calculate ratio between relcov and segment lengths
	# 	#segm$ratio[i] = segm$relcov[i] / (sum(as.numeric(segm$relcov*(segm$end-segm$start))) / sum(as.numeric(segm$end - segm$start)))
	# 	#ix is the segment area from segm applied to the pieces from depcov on the correct chromosome.
	# 	ix <- depcov$chr==segm$chr[i] & (depcov$begin>=segm$start[i] & depcov$end<=segm$end[i])
	# 	segm$ratio[i] <- median(depcov$ratio[ix])
		
	# 	#Index of correct snp positions to pull snp data from
	# 	#                right chromosome
	# 	ix = mastervar$chr==as.character(segm$chr[i]) &
	# 	#    right area (between start and end)
	# 		(mastervar$begin > segm$start[i] & mastervar$begin < segm$end[i])
	# 	segm$snvs[i] = sum(ix)
	# 	segm$ai[i] = 1 - sum(mastervar$min[ix]) / sum(mastervar$max[ix])
	# 	}
	# segm$ai[is.nan(segm$ai)]=0

	#depcov$ai <- NA

	# #En snabbare ai hämtare för depcov segmentstorlek
	# tmp_mv = mastervar[,c(1,2,11,12)]
	# i=1
	# for (j in chroms)
	# 	{
	# 		mv_ix = tmp_mv == j
	# 		mv_ = tmp_mv[mv_ix,]

	# 		while(depcov$chr[i] == j)
	# 			{
	# 				ix = (mv_[,2] >= depcov[i,2] & mv_[,2] <= depcov[i,3])
	# 				depcov[i,6] = 1 - sum(mv_[ix,4]) / sum(mv_[ix,3])
	# 				i = i + 1
	# 				cat(i,"\n")
	# 			}
	# 	}

	# depcov$ai[is.nan(depcov$ai)]=0

#en segare ai hämtare för depcov segmentstorlek (typ en timme)
# for (i in 1:nrow(depcov))
# 	{
# 		ix = mastervar[,1]==depcov[i,1] & (mastervar[,2] >= depcov[i,2] & mastervar[,2] <= depcov[i,3])
# 		depcov[i,6] = 1 - sum(mastervar[ix,12]) / sum(mastervar[ix,11]) #1-min/max
# 		cat(i ,"\n")
# 	}

#Assign AI to depcov segments using min/max from mastervar. (made by Markus)

assignAI <- function(snp_chr,snp_pos,snp_min,snp_max,seg_chr,seg_start,seg_end) {
    nchr <- length(chrs <- unique(seg_chr))
    ai <- NULL
    if (nchr>1) {
        for (c in 1:nchr) {
            ix <- seg_chr==chrs[c]
            ai <- rbind(ai,cbind(chrs[c],assignAI(snp_chr,snp_pos,snp_min,snp_max,seg_chr[ix],seg_start[ix],seg_end[ix])))
        }
        colnames(ai) <- c('chr','begin','end','ai')
        return(ai)
    } else {
        ix <- snp_chr==seg_chr[1]
        snp_pos <- snp_pos[ix]
        snp_min <- snp_min[ix]
        snp_max <- snp_max[ix]
        
        matrix <- rbind(cbind(snp_pos,snp_min,snp_max),cbind(seg_start,-1,-1),cbind(seg_end,-2,-2))
        matrix <- matrix[order(matrix[,2]),]
        matrix <- matrix[order(matrix[,1]),]
        startix <- which(matrix[,2]==-1)
        endix <- which(matrix[,2]==-2)
        
        pos <- matrix[,1]; min <- matrix[,2]; max <- matrix[,3]
        nseg <- length(startix)
        ai=rep(NA,nseg)
        for (i in 1:nseg) {
            if (endix[i]-startix[i] > 1) {
                ai[i] <- 1 - sum(min[ (startix[i]+1):(endix[i]-1) ]) / sum(max[ (startix[i]+1):(endix[i]-1) ])
                #if (ai[i]>1) {
                #    cat(startix[i],endix[i],min[ (startix[i]+1):(endix[i]-1) ],max[ (startix[i]+1):(endix[i]-1) ],'\n')
                #}
            }
        }
        return(data.frame(begin=sort(seg_start),end=sort(seg_end),ai=ai))
    }
}

object = assignAI(mastervar$chr,mastervar$begin,mastervar$min,mastervar$max,
	depcov$chr,depcov$begin,depcov$end)

depcov = merge(depcov,object,by=c("chr","begin","end"))


segm$ai <- segm$snvs <- segm$ratio <- NA

	#Generate the allelic imbalance calculation, among other things.
	for (i in 1:nrow(segm))
		{
		#Calculate ratio between relcov and segment lengths
		#segm$ratio[i] = segm$relcov[i] / (sum(as.numeric(segm$relcov*(segm$end-segm$start))) / sum(as.numeric(segm$end - segm$start)))
		#ix is the segment area from segm applied to the pieces from depcov on the correct chromosome.
		ix <- depcov$chr==segm$chr[i] & (depcov$begin>=segm$start[i] & depcov$end<=segm$end[i])
		segm$ratio[i] <- median(depcov$ratio[ix])
		
			#Index of correct snp positions to pull snp data from
			#                right chromosome
		#ix = mastervar$chr==as.character(segm$chr[i]) &
			#    right area (between start and end)
		#(mastervar$begin > segm$start[i] & mastervar$begin < segm$end[i])
		segm$snvs[i] = sum(ix)
		#segm$ai[i] = 1 - sum(mastervar$min[ix]) / sum(mastervar$max[ix])
		segm$ai[i] <- median(depcov$ai[ix],na.rm=T)
		}
	
	#DNAcopy segmentation of AI

	# library(DNAcopy)
	# data(ideogram,package="patchworkCG")
	# chroms = as.character(unique(ideogram$chr))
	# #chr_list=unique(mastervar$chr)

	# #DNAcopy segmentation of object$ai values.
	# aisegs=NULL
	# for (c in chroms)
	# 	{
	# 	psegments <- qsegments <- NULL
	# 	ix = object$chr==c & object$start < ideogram$start[ideogram$chr==c]
	# 	if (sum(ix)>0)
	# 		{
	# 		psegments = segment(CNA(object$ai[ix], c, object$start[ix],
	# 		 data.type='logratio',sampleid=paste(c,'p')),
	# 		  undo.splits='sdundo', undo.SD=1,min.width=2, alpha=0.0001)$output[,2:6]
	# 		psegments$arm='p'
	# 		}
	# 	ix = object$chr==c & object$start > ideogram$end[ideogram$chr==c]
	# 	if (sum(ix)>0)
	# 		{
	# 		qsegments = segment(smooth.CNA(CNA(object$ai[ix], c, object$start[ix],
	# 		 data.type='logratio',sampleid=paste(c,'q')),smooth.region=40),
	# 		  undo.splits='sdundo', undo.SD=1,min.width=2, alpha=0.0001)$output[,2:6]
	# 		qsegments$arm='q'
	# 		}
	# 	aisegs=rbind(aisegs,psegments,qsegments)
	# 	}
	# colnames(aisegs)=c('chr','start','end','np','mean','arm')


	# 	aisegs$ai <- aisegs$snvs <- aisegs$ratio <- NA

	# 	for (i in 1:nrow(aisegs))
	# 	{
	# 	#Calculate ratio between relcov and aisegsent lengths
	# 	#aisegs$ratio[i] = aisegs$relcov[i] / (sum(as.numeric(aisegs$relcov*(aisegs$end-aisegs$start))) / sum(as.numeric(aisegs$end - aisegs$start)))
	# 	#ix is the aisegsent area from aisegs applied to the pieces from depcov on the correct chromosome.
	# 	ix <- depcov$chr==aisegs$chr[i] & (depcov$begin>=aisegs$start[i] & depcov$end<=aisegs$end[i])
	# 	aisegs$ratio[i] <- median(depcov$ratio[ix])
		
	# 	#Index of correct snp positions to pull snp data from
	# 	#                right chromosome
	# 	ix = mastervar$chr==as.character(aisegs$chr[i]) &
	# 	#    right area (between start and end)
	# 		(mastervar$begin > aisegs$start[i] & mastervar$begin < aisegs$end[i])
	# 	aisegs$snvs[i] = sum(ix)
	# 	aisegs$ai[i] = 1 - sum(mastervar$min[ix]) / sum(mastervar$max[ix])
	# 	}
	
	# #aisegs$ai[is.nan(aisegs$ai)]=0

	# #Convert factor to character
	# aisegs$chr = as.character(aisegs$chr)

	# 	CG_KaCh(aisegs$chr,aisegs$start,aisegs$end,aisegs$ratio,aisegs$ai,
	# 				depcov$chr,depcov$begin,depcov$ratio,(1 - mastervar$min/mastervar$max),
	# 				mastervar$chr,mastervar$begin,
	# 				name=name,xlim=c(0,2.4),ylim=c(0,1))

	cat("Calculations complete \n")
	
	cat("Saving objects to CG.Rdata \n")
	
	save(segm,mastervar,depcov,file="CG.Rdata")
	
	#------------------------------------------------#
	# Plot
	#------------------------------------------------#
	#patchwork.CG.copynumbers(cn2=0.75,delta=0.36,het=0.26,hom=0.97) #HCC1187 100%
	#patchwork.CG.copynumbers(cn2=0.8,delta=0.29,het=0.21,hom=0.8) #HCC1187 70%
	#patchwork.CG.copynumbers(cn2=,delta=,het=,hom=) #HCC2218 100%
	
	cat("Initiating Plotting \n")
	CG_karyotype(segm$chr,segm$start,segm$end,segm$ratio,segm$ai,
					name=name,xlim=c(0,2.4),ylim=c(0,1))
	
	CG_KaCh(segm$chr,segm$start,segm$end,segm$ratio,segm$ai,
					depcov$chr,depcov$begin,depcov$ratio,(1 - mastervar$min/mastervar$max),
					mastervar$chr,mastervar$begin,
					name=name,xlim=c(0,2.4),ylim=c(0,1))
	cat("Plotting Complete \n")
	cat("patchwork.CG.plot Complete","\n","Please read documentation on running patchwork.CG.copynumbers \n")
	}