patchwork.Medians_n_AI <- function(segs,kbsegs,alf)	
	{
	## Get medians and AI for them segs
	segs$median <- segs$ai <- segs$snvs <- NA
	for (i in 1:nrow(segs))
		{
		ix = kbsegs$chr==as.character(segs$chr[i]) & ( kbsegs$pos > segs$start[i] & kbsegs$pos < segs$end[i] )
		segs$median[i] = median(kbsegs$ratio[ix])

		ix2 = alf$achr==as.character(segs$chr[i]) & ( alf$apos > segs$start[i] & alf$apos < segs$end[i] )
		segs$snvs[i] = sum(ix2)
		segs$ai[i] = 1 - sum(alf$amin[ix2]) / sum(alf$amax[ix2])
		#This used for assignAI / updated AI calculation
		#segs$ai[i]=median(kbsegs$ai[ix],na.rm=T)
		}
	segs$ai[is.nan(segs$ai)]=0
	segs = segs[!is.na(segs$ai),]
	return(segs)
	}