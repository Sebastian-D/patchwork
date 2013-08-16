patchwork.alleledata <- function(Pileup, normalalf=NULL, vcf)
	{
	
	#Load data included in package
	packagepath = system.file(package="patchwork")
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))

	#data(commonSnps132,package="patchworkData")

	#if .perl folder exists, delete it.
	#test = getwd()
	#try( setwd(".perl"),silent=TRUE)
	#if(getwd() == paste(test,"/.perl",sep=""))
	#	{
	#	setwd("..")
	#	system("rm -r .perl")
	#	}

	#Copy perl information from package install location
	system(paste("cp -r ",packagepath,"/perl .perl",sep=""))

	if(is.null(vcf))
		{
		#old samtools pileup -vcf used.
		#cat("pileup used \n")
		system(paste("cat ",Pileup," | perl .perl/pile2alleles.pl > ",getwd(),"/pile.alleles",sep=""))
		} else 
		{
		#mpileup and bcftools used. read the pileup and vcf.
		#cat("mpileup used \n")
		system(paste("perl .perl/mpile2alleles.pl ",Pileup," ",vcf," >",getwd(),"/pile.alleles",sep=""))
		}

		#Cleanup
		system("rm -r .perl")

		#Read the generated pile.alleles
		alf = read.csv('pile.alleles', sep='\t',header=F)[,1:6]
		#system("rm pile.alleles")
	
	colnames(alf) = c('achr','apos','atype','aqual','atot','amut')
	alf$aref = alf$atot - alf$amut
	#alf = alf [alf$amut >= 1 & alf$aref >= 1,]
	alf$amax = apply(alf[,6:7],1,max)
	alf$amin = alf$atot-alf$amax
	#alf=alf[alf$atot<100,] # outlierz

	
	#Force compatability between naming of chromosomes to our chromosome names
	#from ideogram file. chr1...chr22,chrX,chrY
	#data(ideogram,package="patchworkData")
	if(length(grep("M",alf$achr))!=0)
		{
		alf = alf[-grep("M",alf$achr),]
		}
	alf$achr = as.character(alf$achr)
	
	x_x = strsplit(alf$achr[1],"chr")
	if(length(x_x[[1]]) == 1)
		{
		for(i in 1:22)
			{
			alf$achr[alf$achr == as.character(i)] = paste("chr",as.character(i),sep="")
			}
		i = "X"
		alf$achr[alf$achr == as.character(i)] = paste("chr",as.character(i),sep="")
		i = "Y"
		alf$achr[alf$achr == as.character(i)] = paste("chr",as.character(i),sep="")
		}
	
	#Check if HG18 or HG19 should be applied
	con = file(Pileup,"rt")
	line = readLines(con,1)
	hgcheck = strsplit(line,"\t")

	#Sometimes chr1 isnt first in the pileup file, traverse until you get to chr1 to see if hg18 or hg19 should be applied.
		#while hgcheck is not "chr1" AND not "1"
	while((hgcheck[[1]][1] != "chr1") && (hgcheck[[1]][1] != "1"))
		{
		line = readLines(con,1)
		hgcheck = strsplit(line,"\t")
		}
	close(con)

	#Check pos to see which Snps should be used, hg18 or hg19
	if (as.numeric(hgcheck[[1]][2]) <= 10000)
		{
		data(commonSnpsHG18,package="patchworkData")
  } else
		{
		data(commonSnps132,package="patchworkData")
		}

	# Annotate with dbSNP
	dbSnp=dbSnp[,c(1,3)]; dbSnp$dbSnp=T
	alf=merge(alf,dbSnp, all.x=T, all.y=F, by=1:2)
	alf$dbSnp[is.na(alf$dbSnp)]=F

	#Finally, if a matched normal was available, remove SNPs that were not (somewhat) heterozygous there.
	#if (!is.null(normalalf)) {
	#	normalalf <- normalalf[normalalf$amin/normalalf$atot > 0.1,]
	#	alf <- merge(normalalf[,1:2],alf,by=1:2,all=F)
	#}
	
	return(alf)
	}
