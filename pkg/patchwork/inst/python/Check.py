import sys
import os

try:
    import pysam
    print "X"
except:
	#If import pysam failed
    print "Pysam Not Available"
    sys.exit()
