import sys                                                                       
import string
import re                                                                        
                                                                                 
#This program  reads in the dot file produced by valgrind 
# and it ignores some of library nodes (e.g. libblas)

f = open(sys.argv[1], 'rU')               
lines = f.readlines()

fileStr = []
for ln in range(0, len(lines)):  ## iterates over the lines of the file
    line = lines[ln]
    a = string.split(line)
    if (string.find(a[0], "\"") != -1):
	continue
    else:
	if (len(a)>1 and string.find(a[1], "->") != -1):
	  if (a[0][0] < 'A' or a[0][0] > 'Z') or (a[2][0] < 'A' or a[2][0] > 'Z'):
	    continue
	print line,

