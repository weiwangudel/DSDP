import sys                                                                       
import string
import re                                                                        
                                                                                 
#This program  reads in the dot file produced by valgrind 
# and it ignores some of library nodes (e.g. libblas)


# read in the csdp function names
f_ds = open(sys.argv[2], 'r')
all_csdp = f_ds.readlines()


def belongs_to(name, names):
  for i in range (0, len(all_csdp)):
    if (string.find(names[i], name) != -1):
      return True
  return False

# read in the big dot to filter
f = open(sys.argv[1], 'rU')               
lines = f.readlines()

fileStr = []
for ln in range(0, len(lines)):  ## iterates over the lines of the file
    line = lines[ln]
    a = string.split(line)
    if (len(a)>1 and string.find(a[1], "->") != -1):
	if (not belongs_to (a[0], all_csdp)) or (not belongs_to (a[2], all_csdp)):  # ignore non-csdp edges
	    continue
	print line,
    elif ((string.find(line, "->") == -1) and (string.find(line, sys.argv[3]) == -1) and (string.find(line, "label") != -1))  :     # ignore non-csdp nodes 
	continue
    else:
	print line,

