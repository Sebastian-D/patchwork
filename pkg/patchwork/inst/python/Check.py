import sys
import os

try:
    import pysam
    print "X"
    sys.exit()
except:
	#It didnt work.
    print "Pysam Not Available"
    sys.exit()
