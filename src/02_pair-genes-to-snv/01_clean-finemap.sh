#!/bin/bash

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

# Output files: 
# fm_nc.txt: list of non-converning loci, that is loci where cumulative post. prob. of k from 1-19 is below 95%.
# fm_mhc.bed: list of loci that overlap the MHC region on hg19.
# fm_rem.bed: a non-redundant list of loci that overlap the MHC or failed to converge.
# fm_f1.txt: List of all FINEMAP results, with locus id appended as first column.
# fm_f2.txt: List of all FINEMAP results from fm_f1.txt, excluding SNPs that overlap regions flagged for removal (non-converged and MHC).
# fm_f3.txt: List of all FINEMAP results from fm_f2.txt, excluding SNPs with post. prob == 0. 

set -eux
# Specify directory containing FINEMAP analysis results. 
# Likely ../finemap as the current directory is typically in same dir as the FINEMAP results.
# FINEMAP_DIR=../finemap
FINEMAP_DIR=${1:-../finemap}
# FINEMAP_DIR=$1

# List non-converged loci, ie those with cumulative posterior probability < 0.95 for range of 1 to 19 SNPs being causal. 
# Not that when FINEMAP fails the *.log_sss file is empty or incomplete. 
# It this case I initialized the cumulative posterior probability to be 0, thus these loci will also be removed.

for f in ${FINEMAP_DIR}/regions/L*.log_sss
do
    locus=$(basename $f .log_sss)
    grep -A 20  Post-Pr $f | tail -n +2 | perl -p -e "s/.*-> ([0-9.e-]+)\)?/\1/" \
	| awk -vlocus=$locus 'BEGIN { s = 0 } { s += $1 } END { if (s < 0.95 ) { print locus, s }}'

done > fm_nc.txt

# Some old code that find where post-prob at k >= 19 is maximum prob and removes these loci.
# grep -A 21  Post-Pr $f | tail -n +2 | perl -p -e "s/.*-> ([0-9.e-]+)\)?/\1/" | awk -vlocus=$locus 'BEGIN { max=0 } $1 > max { max=$1; k=NR-1 } END { if (k >= 19) { print locus, k, max } }'

# Find loci that overlap the MHC locus
bedtools intersect -b ../../mhc.bed -a ${FINEMAP_DIR}/leadLoci.bed > fm_mhc.bed

# Get listing of all loci to remove, ie those that failed FINEMAP convergence and overlap the MHC locus.
(cut -f 1 -d ' ' fm_nc.txt; cut -f 4 fm_mhc.bed) | sort | uniq > fm_rem.txt

# Remove loci from the loci bed file
awk 'BEGINFILE { FC += 1 } FC==1 { m[$1]=$1;next} !($4 in m) { print }' fm_rem.txt ${FINEMAP_DIR}/leadLoci.bed > leadLoci.bed
# awk 'FNR==NR { m[$1]=$1;next} !($4 in m) { print }' fm_rem.txt ${FINEMAP_DIR}/leadLoci.bed > leadLoci.bed
# Merge overlapping loci
bedtools merge -i leadLoci.bed -c 4 -o collapse -delim '|' | sort -k1,1 -k2,2n > leadLoci_merged.bed

# Retain lead SNPs from remaining loci
bedtools intersect -a leadLoci.bed -b ${FINEMAP_DIR}/leadSnps.tab -wa -wb | cut -f 5-8 | sort -k1,1 -k2,2n | uniq > leadSnps.bed

# Merge all FINEMAP results, appending locus name as first column
awk 'FNR==1 { match(FILENAME, /(L[0-9]+)/,a); locus=a[1] } FNR==1 && NR == 1 { print "locus", $0; next } FNR==1 { next } { print locus, $0 }' \
    ${FINEMAP_DIR}/regions/L*.snp \
    > fm_f1.txt

# Remove loci flagged for removal
awk 'BEGINFILE { FC += 1 } FC==1{ removals[$1]=$1; next } !($1 in removals) { print }' \
    fm_rem.txt \
    fm_f1.txt \
    > fm_f2.txt

# Retain lead SNPs from FINEMAP results, retaining lowest log10bf where SNP is present more than once due to loci overlap.
(head -n 1 fm_f2.txt; awk 'FNR==NR { m[$4]=$4; next } $3 in m { print }' leadSnps.bed fm_f2.txt | sort -k 3,3 -k13,13n | sort -k3,3 -u | sort -k4,4 -k5,5n) > fm_leadSnps.txt

# Remove SNPs with FINEMAP post. prob. < 0, which is equivalent to a log10BF of "-inf".
awk '$12 > 0 { print }' fm_f2.txt > fm_f3.txt

# Count duplicate SNPs, report tally of times SNP was present
tail -n +2 fm_f3.txt | cut -f 3 -d ' ' | sort | uniq -c > fm_tally.txt

echo " loci:" $(awk 'END { print NR }' ${FINEMAP_DIR}/leadLoci.bed) $(ls ${FINEMAP_DIR}/regions/L*.log_sss | awk 'END { print NR }') $(tail -n +2 fm_f1.txt | awk '{ m[$1]=$1 } END { print length(m) }') > 00_clean-finemap.log
echo "  rem:" $(awk 'END { print NR }' fm_nc.txt) $(awk 'END { print NR }' fm_mhc.bed) $(awk 'END { print NR }' fm_rem.txt) >> 00_clean-finemap.log
echo "final:" $(awk 'END { print NR }' leadLoci.bed) $(tail -n +2 fm_f2.txt | awk '{ m[$1]=$1 } END { print length(m) }') $(awk 'END { print NR }' leadLoci_merged.bed) >> 00_clean-finemap.log
echo " lsnp:" $(awk 'END  { print NR }' ${FINEMAP_DIR}/leadSnps.bed) $(awk 'END { print NR }' leadSnps.bed) $(awk 'END { print NR-1 }' fm_leadSnps.txt) >> 00_clean-finemap.log
echo "fmsnp:" $(awk 'END { print NR }' fm_f1.txt) $(awk 'END { print NR }' fm_f2.txt) $(awk 'END { print NR }' fm_f3.txt) >> 00_clean-finemap.log
echo "tally:" $(awk '{ print $1 }' fm_tally.txt | sort | uniq -c | sort -k2,2n | tr '\n' ' ') >> 00_clean-finemap.log


# Sort remaining results by SNP then ascending log10bf, retain first result if SNP duplicated (lowest log10bf).

(head -n 1 fm_f3.txt; tail -n +2 fm_f3.txt | sort -k 3,3 -k13,13n | sort -k3,3 -u | sort -k4,4 -k5,5n) > fm.txt
awk '$13 >= 1' fm.txt > fm1.txt
awk '$13 >= 2' fm.txt > fm2.txt
awk '$13 >= 3' fm.txt > fm3.txt

# Convert to BED format. Note position is 0-based.
# locus index rsid             chromosome position  allele1 allele2 maf      beta        se         z         prob        log10bf   group corr_group prob_group  log10bf_group mean         sd          mean_incl   sd_incl
# L635  1160  10:101769738_T_G 10         101769738 T       G       0.007825 -0.0490727  0.0138316  -3.54787  0.00826392  0.879261  1160  1          0.00826392  0.879261      -0.000363004 0.00412466  -0.0439264  0.0120459

awk '{ n=$1","$3; printf "chr%d\t%d\t%d\t%s\t%s\t%s\n", $4, $5-1, $5, n, $13, "+" }' <(tail -n+2 fm.txt)  | sort -k1,1 -k2,2n > fm.bed
awk '{ n=$1","$3; printf "chr%d\t%d\t%d\t%s\t%s\t%s\n", $4, $5-1, $5, n, $13, "+" }' <(tail -n+2 fm1.txt) | sort -k1,1 -k2,2n > fm1.bed
awk '{ n=$1","$3; printf "chr%d\t%d\t%d\t%s\t%s\t%s\n", $4, $5-1, $5, n, $13, "+" }' <(tail -n+2 fm2.txt) | sort -k1,1 -k2,2n > fm2.bed
awk '{ n=$1","$3; printf "chr%d\t%d\t%d\t%s\t%s\t%s\n", $4, $5-1, $5, n, $13, "+" }' <(tail -n+2 fm3.txt) | sort -k1,1 -k2,2n > fm3.bed
awk '{ n=$1","$3; printf "chr%d\t%d\t%d\t%s\t%s\t%s\n", $4, $5-1, $5, n, $13, "+" }' <(tail -n+2 fm_leadSnps.txt) | sort -k1,1 -k2,2n > fml.bed

echo " fmbf:" 0:$(awk 'END { print NR }' fm.bed) 1:$(awk 'END { print NR }' fm1.bed) 2:$(awk 'END { print NR }' fm2.bed) 3:$(awk 'END { print NR }' fm3.bed) >> 00_clean-finemap.log



