#!/bin/bash

python ~/Desktop/gprof2dot.py  -f callgrind -e 5 -n 5 -z DSDPSolve $1  > $1.dot
