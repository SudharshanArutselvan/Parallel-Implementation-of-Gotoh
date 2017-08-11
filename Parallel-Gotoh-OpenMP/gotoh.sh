#!/bin/bash
#$ -N Gotoh
#$ -q class64-amd
#$ -pe openmp 8

# Module load gcc compiler version 4.9.2
module load  gcc/4.9.2

# Runs a bunch of standard command-line
# utilities, just as an example:

echo "Script began:" `date`
echo "Node:" `hostname`
echo "Current directory: ${PWD}"

echo ""
for trial in 1 2 3 4 5; do
  echo "*** Trial ${trial} ***"
	./align sequence1.fasta sequence2.fasta
done

echo ""
echo "=== Done! ==="

# eof
