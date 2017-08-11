#!/bin/bash
#$ -N ADD
#$ -q gpu2

# Notes: 
# The GPU node (compute-1-14) has 4 Tesla M2090 GPU cards. The node also
# has 24 Intel cores.
#
# In the Grid Engine directives above, we are requesting one (1) GPU card:
#    -l gpu=1
# and 6 Intel cores with the Parallel Environment "gpu-node-cores":
#    -pe gpu-node-cores 6   
#
# When requesting GPU resources, please try requesting 6 Intel cores per
# each gpu card you request.  Since the node has 24 Intel cores, the
# division comes out to 6 Intel cores per each GPU card.  
#
# There are no fixed numbers when requesting cores verses GPU cards, it all
# depends on the running program. If you can run with 2 Intel cores and 2 
# GPU cards, then use those numbers.
#
# To see how many GPU cards compute-1-14 has avaialble for a job:
# qhost -F gpu -h compute-1-14
#
# For details, see: http://hpc.oit.uci.edu/gpu
#

# Module load Cuda Compilers and GCC
module load  cuda/8.0
module load  gcc/4.4.3

# Runs a bunch of standard command-line
# utilities, just as an example:

echo "Script began:" `date`
echo "Node:" `hostname`
echo "Current directory: ${PWD}"

nvcc gotoh.cc timer.c  main5.cu -o main
./main sequence3.fasta sequence4.fasta -gapopen 2 -gapextend 1

echo ""
echo "=== Done! ==="

# eof
