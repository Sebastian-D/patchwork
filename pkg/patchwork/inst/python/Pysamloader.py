import sys
import re
import os

import pysam

try:
    infile = pysam.Samfile(sys.argv[1],"rb")
except IOError,eStr:
    print "Error:Cannot open",infile," for reading: ", eStr
    sys.exit()
    
outfile = open(".Tmp_chr_pos","w")

try:
	for read in infile.fetch(sys.argv[2]):
		outfile.write(str(read.pos) + "\n")
except:
	a = re.split('chr',sys.argv[2])
	for read in infile.fetch(a[1]):
		outfile.write(str(read.pos) + "\n")
	
infile.close()
outfile.close()
