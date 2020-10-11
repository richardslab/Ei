#!/bin/bash

## GWAS summary statistics file
GWASF=gwas.sumstats.txt

## Example index and name of columns from GWAS summary statistics file (useful to map to output columns in AWK code below). This will be specific to each input GWAS.
##    1       2       3      4       5      6         7        8     9     10
##  CHR     SNP     POS     A1      A2      N       AF1     beta    se      p


## Header of 2 output files, one used for PLINK LD clumping and the other for FINEMAP.
echo "SNP CHR BP A1 A2 FRQ BETA SE P" > gwas.assoc
echo "rsid chromosome position allele1 allele2 eaf beta se p" > gwas.txt

## AWK code harmonizes GWAS sumstats to include only SNPs with MAF > 0.005 and those present in UKB, matched by position and alleles. Not that the SNP id used is the rs# provided by UKB. The format of the UKB SNP index is provided below.
awk 'FNR==NR { m[$2]=$1; next } $7 > 0.995 || $7 < 0.005 { next } ($1":"$3":"$4":"$5 in m) { OFS=" "; print m[$1":"$3":"$4":"$5],$1,$3,$4,$5,$7,$8,$9,$11; next } ($1":"$3":"$5":"$4 in m) { OFS=" "; print m[$1":"$3":"$5":"$4],$1,$3,$4,$5,$7,$8,$9,$11 }' \
    ../ukb.snps.txt \
    ${GWASF} \
    | tee -a gwas.txt >> gwas.assoc

## $ head ../ukb.snps.txt
## rs367896724 1:10177:A:AC
## rs540431307 1:10235:T:TA
## rs201106462 1:10352:T:TA
## rs548419688 1:10505:A:T
## rs568405545 1:10506:C:G
## rs534229142 1:10511:G:A