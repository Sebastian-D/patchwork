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

assignAI <- function(snp_chr,snp_pos,snp_min,snp_max,seg_chr,seg_start,seg_end) 
{
    nchr <- length(chrs <- unique(seg_chr))
    ai <- NULL
    if (nchr>1) 
      {
        for (c in 1:nchr) 
          {
            ix <- seg_chr==chrs[c]
            ai <- rbind(ai,cbind(chrs[c],assignAI(snp_chr,snp_pos,snp_min,snp_max,seg_chr[ix],seg_start[ix],seg_end[ix])))
          }
        colnames(ai) <- c('chr','pos','end','ai')
        return(ai)
      } 
      else 
        {
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
        for (i in 1:nseg) 
          {
            if (endix[i]-startix[i] > 1) 
              {
                ai[i] <- 1 - sum(min[ (startix[i]+1):(endix[i]-1) ]) / sum(max[ (startix[i]+1):(endix[i]-1) ])
                #if (ai[i]>1) {
                #    cat(startix[i],endix[i],min[ (startix[i]+1):(endix[i]-1) ],max[ (startix[i]+1):(endix[i]-1) ],'\n')
                #}
              }
          }
        return(data.frame(begin=sort(seg_start),end=sort(seg_end),ai=ai))
        }
}

###################################