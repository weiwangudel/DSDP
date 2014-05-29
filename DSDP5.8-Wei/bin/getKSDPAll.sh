#!/bin/bash

python gprof2dot.py  -f callgrind -e 0 -n 0 -z KSDPConeComputeHessian $1  > $1.dot
