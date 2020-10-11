#!/bin/bash

# this script relied on code from Makefile
# Build lead SNPs file from results of LD clumping
make leadSnps.bed
# Using leadSnps.bed build a set of loci
make leadLoci.bed


# the code below should be looped over all loci from leadLoci.bed. Example is provided for locus "L1" from chromosome 1.

LOCUS=L1
CHROM=1
export $LOCUS
export $CHROM
echo $LOCUS $CHROM
# Generate FINEMAP input files
make finemap LOCUS=${LOCUS} CHROM=${CHROM}
cd regions/
# Run FINEMAP using inpute files from previous step
make -f ../Makefile run_finemap LOCUS=${LOCUS} CHROM=${CHROM}
