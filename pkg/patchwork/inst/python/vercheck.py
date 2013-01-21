#!/usr/bin/python

# import libraries
import urllib2
import re
 
  
### MAKE SURE THIS URL IS CORRECT
# fetch the html from the package's r-forge page
html = urllib2.urlopen('https://r-forge.r-project.org/R/?group_id=1250', timeout = 2).read()
   
    
# parse the html for the version number
match = re.search('Version: <b>(.+)</b>', html)
     
# print the version number if found
if match:
        print match.groups()[0]

