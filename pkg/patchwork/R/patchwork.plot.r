patchwork.plot <- function(Tumor.bam,Tumor.pileup,Tumor.vcf=NULL,Normal.bam=NULL,Normal.pileup=NULL,Normal.vcf=NULL,
							Reference=NULL,Alpha=0.0001,SD=1)
	{

	#Extract the name from tumor.bam

	name = strsplit(tolower(Tumor.bam),".bam")

	name = strsplit(name[[1]],"/")

	name = name[[1]][length(name[[1]])]

	alf <- normalalf <- NULL

	#Main check function. If this passes, and name_copynumbers.Rdata exists
	#it will jump straight to plotting.
	try( load(paste(name,"_copynumbers.Rdata",sep="")),silent=TRUE)

	if(is.null(alf))
		{
		#Attempt to read file pile.alleles, incase its already been created.
		try( load(paste(name,"_pile.alleles.Rdata",sep="")), silent=TRUE )
	
		#If it wasnt created, create it using the perl script pile2alleles.pl
		#which is included in the package. Creates pile.alleles in whichever folder
		#you are running R from. (getwd())
		if(is.null(alf))
			{
			cat("Initiating Allele Data Generation \n")
			if (!is.null(Normal.pileup)) normalalf <- patchwork.alleledata(Normal.pileup, vcf=Normal.vcf)
			alf = patchwork.alleledata(Tumor.pileup, normalalf=normalalf,Tumor.vcf)
			cat("Allele Data Generation Complete \n")
			save(alf, normalalf,file=paste(name,"_pile.alleles.Rdata",sep=""))
			}
	
		if (!is.null(normalalf)) 
			{
			normalalf <- normalalf[normalalf$amin/normalalf$atot > 0.2,]
			alf <- merge(normalalf[,1:2],alf,by=1:2,all=F)
			}




		#Generate chroms object depending on the naming in pileup (alf).
		#either chr1...chr22,chrX,chrY or 1...22,X,Y
		data(ideogram,package="patchworkData")
	
		chroms = as.character(unique(ideogram$chr))	
	
		## Read coverage data. If already done, will load "data.Rdata" instead.
		data <- normaldata <- NULL
		try ( load(paste(name,"_data.Rdata",sep="")), silent=TRUE )
	
		#If data object does not exist, IE it was not loaded in the previous line
		#perform the function on the chromosomes to create the object.
		if(is.null(data)) 
			{
			cat("Initiating Read Chromosomal Coverage \n")
			data = patchwork.readChroms(Tumor.bam,chroms)
			cat("Read Chromosomal Coverage Complete \n")
			save(data,file=paste(name,"_data.Rdata",sep=""))
			}
	
		#Perform GC normalization if the amount of columns indicate that gc normalization
		#has not already been performed.
		if(length(data) < 7) 
			{
			cat("Initiating GC Content Normalization \n")
			data = patchwork.GCNorm(data)
			cat("GC Content Normalization Complete \n")
			save(data,file=paste(name,"_data.Rdata",sep=""))
			}
		#after this, no further normalization was done on reference sequences.

		#If matched normal bam file was defined, load and normalize it the same way.
		if (!is.null(Normal.bam) & is.null(Reference)) 
			{
			normaldata=NULL
			try ( load(paste(name,"_normaldata.Rdata",sep="")), silent=TRUE )
			if(is.null(normaldata)) 
				{
				cat("Initiating Read Chromosomal Coverage (matched normal) \n")
				normaldata = patchwork.readChroms(Normal.bam,chroms)
				cat("Read Chromosomal Coverage Complete (matched normal) \n")
				}	
			if(length(normaldata) != 3) 
				{
				cat("Initiating GC Content Normalization (matched normal) \n")
				normaldata = patchwork.GCNorm(normaldata[,1:6])
				cat("GC Content Normalization Complete (matched normal) \n")
				normaldata <- data.frame(chr=normaldata$chr,pos=normaldata$pos,normal=normaldata$norm)
				save(normaldata,file=paste(name,"_normaldata.Rdata",sep=""))
				}
			#save(normaldata,file='normaldata.Rdata')
			} 

	
		# Smooth the data.
		kbsegs = NULL
		try( load(paste(name,"_smoothed.Rdata",sep="")), silent=TRUE )

		if(length(kbsegs) == 0) 
			{
			cat("Initiating Smoothing \n")
			kbsegs = patchwork.smoothing(data,normaldata,Reference,chroms)
			save(kbsegs,file=paste(name,"_smoothed.Rdata",sep=""))
			cat("Smoothing Complete \n")
			}
	
		#Segment the data.
		library(DNAcopy)
		segs = NULL
		try( load(paste(name,"_Segments.Rdata",sep="")), silent=TRUE )
		if(length(segs) == 0)
			{
			cat("Initiating Segmentation \n")
			cat("Note: If segmentation fails to initiate the probable reason is that you have not ")
			cat("installed the R package DNAcopy. See the homepage, http://patchwork.r-forge.r-project.org/ ,
				or ?patchwork.readme for installation instructions. \n")
			segs = patchwork.segment(kbsegs,chroms,Alpha,SD)
			save(segs,file=paste(name,"_Segments.Rdata",sep=""))
			cat("Segmentation Complete \n")
			}

		#Assign AI to kbsegs. (not in use atm)   ---MARKUS ---
		#object = assignAI(alf$achr,alf$apos,alf$amin,alf$amax,
		#kbsegs$chr,(kbsegs$pos-5000),(kbsegs$pos+5000))

		#kbsegs = cbind(kbsegs[,1:5],object[,4])
		#colnames(kbsegs)=c("chr","pos","coverage","refcoverage","ratio","ai")
		
		if(length(segs) == 6)
			{
			#Get medians and AI for the segments.
			cat("Initiating Segment data extraction (Medians and AI) \n")
			segs = patchwork.Medians_n_AI(segs,kbsegs,alf)
			save(segs,file=paste(name,"_Segments.Rdata",sep=""))
			cat("Segment data extraction Complete \n \n \n")
			}
	
		cat(paste("Saving information objects needed for patchwork.copynumbers() in ",name,"_copynumbers.Rdata \n \n \n",sep=""))
		save(segs,alf,kbsegs,file=paste(name,"_copynumbers.Rdata",sep=""))
		}
	
	#Plot it
	cat("Initiating Plotting \n")
	karyotype(segs$chr,segs$start,segs$end,segs$median,segs$ai,
			name=as.character(name),
			xlim=c(0,2.4),
			ylim=c(0.1,1))
	karyotype_chroms(segs$chr,segs$start,segs$end,segs$median,segs$ai,
			kbsegs$chr,kbsegs$pos,kbsegs$ratio,
			alf$achr,alf$apos,(1-alf$amin/alf$amax),
			name=as.character(name),
			xlim=c(0,2.4),
			ylim=c(0.1,1))
	cat("Plotting Complete \n")
	cat("patchwork.plot Complete.\n")
	cat("Below you may see some warning messages, you can read about these on our homepage. They are nothing to be worried about. \n")
	}
	
	