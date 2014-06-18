#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "Usage:"
  echo $0 "csdp|dsdp5 [threshold:e.g. 0.01]"
  exit 1 
fi

if [ $# -eq 1 ]; then 
  thresh=0             #default threshold
else 
  thresh=$2            #user specified threshold
fi

# 3 inputs: input_100.txt, input_1317.txt, input 9955.txt
for i in 100 1317 9955 
do 
  echo "Running valgrind for input_$i.txt" 
               #echo "simulation result is saved to /dev/null"
               #echo "valgrind output is saved to callgrind.out.$i"

  #valgrind --tool=callgrind --callgrind-out-file=callgrind.out.$i ./$1  input_$i.txt   > /dev/null 2>/dev/null
  echo "Got Callgrind Output, generating dot file for input_$i.txt"
  python gprof2dot.py  -f callgrind -e 0 -n $thresh -z main callgrind.out.$i  > $i.dot

  #get csdp|dsdp5 only functions
  grep $1 $i.dot |cut -d" " -f1 |awk '{print $1}'  > $i-names
  python ignore-some-library.py $i.dot $i-names $1 > $i-$1-only.dot
  #dot -Tpdf $i-$1-only.dot -o $i-$1-only.pdf
  rm $i-names
  rm $i.dot
done 


  echo "Callgraphs (dot) generated, combining into one......."
  # need to replace this part if we want to strictly distingush inputs
  sed -i 's/0d0d73/000000/g'  1317-$1-only.dot    #change color
  sed -i 's/0d0d73/008000/g'  9955-$1-only.dot    #change color
  
  #preparing to concatenate them 
  #get rid of the first 4th and/or the last line
  head -n -1 100-$1-only.dot > three_in_one.dot 
  tail -n +5 1317-$1-only.dot | head -n -1 >> three_in_one.dot 
  tail -n +5 9955-$1-only.dot  >> three_in_one.dot 
  dot -Tpdf three_in_one.dot -o three_in_one.pdf 
  
  echo "Callgraphs three in one generated: blue - black - green "  
