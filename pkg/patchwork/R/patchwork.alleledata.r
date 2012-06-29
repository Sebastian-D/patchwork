patchwork.alleledata <- function(Pileup, normalalf=NULL)
	{
	
	#Load data included in package
	packagepath = system.file(package="patchwork")
	#load(paste(packagepath,"/data/commonSnps132.RData",sep=""))
	
	data(commonSnps132,package="patchworkData")

	#if .perl folder exists, delete it.
	#test = getwd()
	#try( setwd(".perl"),silent=TRUE)
	#if(getwd() == paste(test,"/.perl",sep=""))
	#	{
	#	setwd("..")
	#	system("rm -r .perl")
	#	}
	
	system(paste("cp -r ",packagepath,"/perl .perl",sep=""))
	system(paste("cat ",Pileup," | perl .perl/pile2alleles.pl > ",getwd(),"/pile.alleles",sep=""))
	system("rm -r .perl")
	alf = read.csv('pile.alleles', sep='\t',header=F)[,1:6]
	#system("rm pile.alleles")
	
	colnames(alf) = c('achr','apos','atype','aqual','atot','amut')
	alf$aref = alf$atot - alf$amut
	alf = alf [alf$amut >= 1 & alf$aref >= 1,]
	alf$amax = apply(alf[,6:7],1,max)
	alf$amin = alf$atot-alf$amax
	#alf=alf[alf$atot<100,] # outlierz

	
	#Force compatability between naming of chromosomes to our chromosome names
	#from ideogram file. chr1...chr22,chrX,chrY
	data(ideogram,package="patchworkData")
	alf = alf[-grep("M",alf$achr),]
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
	
	
	alf=merge(alf,dbSnp[,c(1,3)], all=F, by=1:2)

	#Finally, if a matched normal was available, remove SNPs that were not (somewhat) heterozygous there.
	if (!is.null(normalalf)) {
		normalalf <- normalalf[normalalf$amin/normalalf$atot > 0.1,]
		alf <- merge(normalalf[,1:2],alf,by=1:2,all=F)
	}
	
	return(alf)
	}
