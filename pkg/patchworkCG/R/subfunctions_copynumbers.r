###################### functions
weightedMedian <- function(data,weights) {
    ix <- rep(T,length(data))
    if (length(ix)==0) return(NULL)
    try( ix <- !is.na(data+weights) & !is.nan(data+weights), silent=T)
    try( return (median(rep(data[ix],ceiling(weights[ix])))), silent=T)
    return (NULL)
}
weightedMean <- function(data,weights) {
    ix <- rep(T,length(data))
    if (length(ix)==0) return(NULL)
    try( ix <- !is.na(data+weights) & !is.nan(data+weights), silent=T)
    try( return (mean(rep(data[ix],ceiling(weights[ix])))), silent=T)
    return (NULL)
}
is.autosome <- function(vector) {
    temp <- deChrom_ucsc(vector)
    return (temp<23)
}
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

###################################