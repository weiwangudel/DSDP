#!/bin/bash

python ~/Desktop/gprof2dot.py  -f callgrind -e 0 -n 0 -z DSDPSolve $1  > $1.dot
