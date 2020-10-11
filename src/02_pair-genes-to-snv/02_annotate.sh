#!/bin/bash
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

set -euxo pipefail

TRAIT=$1

BF=2
FMP=fm${BF}
FMB=${FMP}.bed

## Genes that overlap each locus
# find overlapping genes, need to overlap > 50% of the gene to avoid genes being duplicated due to "linking" adjecent yet non-overlapping loci

bedtools intersect -a leadLoci_merged.bed -b ../../../../data/gene_sets/gencode29b37_protein_coding.ncbi_gene_id.bed -wa -wb -F 0.5 > loci_gene_overlaps.bed
awk '{ OFS="\t"; print $5,$6,$7,$8","$4,$9,$10 }' loci_gene_overlaps.bed | sort -k1,1 -k2,2n > loci_genes.bed

# TSS of genes that overlap the locus
awk '$6=="+" { tss=$2 } $6=="-" { tss=$3 } { OFS="\t"; print $1,tss,tss,$4,$5,$6 }' loci_genes.bed | sort -k1,1 -k2,2n > loci_tss.bed 

# 1. Annotate FINEMAP SNPs for overlap with gene body.
bedtools intersect -a loci_genes.bed -b ${FMB} -wa -wb > ${FMP}_gene.bed

# 2. Annotate FINEMAP SNPs for those closest to 5' end of gene.
bedtools closest -a loci_tss.bed -b ${FMB} -wa -wb -d > ${FMP}_tss.bed

# 2.1 Annotate FINEMAP SNPs for those closest to 5' end of gene, upstream only.
# bedtools closest -a loci_tss.bed -b ${FMB} -fu -D a -wa -wb > ${FMP}_upstream.bed


# 3. Annotate FINEMAP SNPs for functional consequence in protein coding gene using dbSNP
bedtools intersect -f 1 -r -a /mnt/RICHARDS_JBOD1/SHARE/DATA/UCSC/hg19/database/snp150_func.bed -b ${FMB} -wa -wb > ${FMP}_dbsnp.txt
awk '$9 ~ $4 || $4 !~ /rs/ { OFS="\t"; print $1,$2,$3,$9,$5,"+" }' ${FMP}_dbsnp.txt > ${FMP}_dbsnp.bed
bedtools intersect -a loci_genes.bed -b ${FMP}_dbsnp.bed -wa -wb > ${FMP}_func.bed

# 4. Annotate FINEMAP SNPs for functional consequence in protein coding gene using SNPEff
awk 'BEGIN { OFS="\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER" } { print $4,$5,$1","$3,$6,$7,".","." }' <(tail -n +2 ${FMP}.txt) > ${FMP}.vcf
java -Xmx4g -jar ~/snpEff/snpEff.jar -v GRCh37.75 ${FMP}.vcf > ${FMP}.ann.vcf 2> ${FMP}.ann.log
grep -v '^#' ${FMP}.ann.vcf | awk '{ split($8,a,","); for (i=1;i<=length(a);i+=1) { split(a[i],b,"|"); OFS="\t"; print $3,b[3],b[4],b[5] }}' | sort | uniq >  ${FMP}_gene_snpeff.txt
awk 'BEGIN { FN=0 } FNR==1 { FN++ } FN==1 { snps[$4]=$0 } FN==2 { split($4,a,","); genes[a[3]]=$0 } FN==3 && ($3 in genes){ OFS="\t"; split(snps[$1],a,"\t"); print genes[$3], a[1],a[2],a[3],a[4],$2,a[6] }' ${FMP}.bed loci_genes.bed ${FMP}_gene_snpeff.txt > ${FMP}_snpeff.bed
if [ -s "${FMP}_snpeff.bed" ]
then
    awk '{ OFS="\t"; print $1,$2,$3,$10":"$4,$11 }' ${FMP}_snpeff.bed | sort -k1,1 -k2,2n | bedtools groupby -o distinct -g 1-4 -c 5 > ${FMP}_snpeff_distinct.bed
else
    touch ${FMP}_snpeff.bed 
fi
# 5. Annotate FINEMAP SNPs for overlap with one or more openchromatin peaks
bedtools intersect -b ../../MAURANO-DHS/${TRAIT}_dhs.bed -a ${FMP}.bed -wa -wb > ${FMP}_opchr.txt
awk '{ OFS="\t"; print $1,$2,$3,$4,$10,"+" }' ${FMP}_opchr.txt > ${FMP}_openchromfull.txt
bedtools groupby -i ${FMP}_openchromfull.txt -o count_distinct -g 1-4 -c 5 | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,"+" }' > ${FMP}_openchromgrp.txt
bedtools closest  -a loci_tss.bed -b ${FMP}_openchromgrp.txt -wa -wb > ${FMP}_openchrom.bed
bedtools closest -t first -b loci_tss.bed -a ${FMP}_openchromgrp.txt -wa -wb > ${FMP}_openchromsnp.tmp
paste <(cut -f 7-12 ${FMP}_openchromsnp.tmp) <(cut -f 1-6 ${FMP}_openchromsnp.tmp) > ${FMP}_openchromsnp.bed

ALL_DHS=$SHARE/DATA/FUNCTIONAL_GENOMICS/OPEN_CHROMATIN/ENCODE_REMC_FLER_T2D_2019aug.flt.hotspots2.fdr0.05.hg19.bed
bedtools intersect -b ${ALL_DHS} -a ${FMP}.bed -wa -wb > ${FMP}_allopchr.txt
awk '{ OFS="\t"; print $1,$2,$3,$4,$10,"+" }' ${FMP}_allopchr.txt > ${FMP}_allopenchromfull.txt
bedtools groupby -i ${FMP}_allopenchromfull.txt -o count_distinct -g 1-4 -c 5 | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,"+" }' > ${FMP}_allopenchromgrp.txt
bedtools closest  -a loci_tss.bed -b ${FMP}_allopenchromgrp.txt -wa -wb > ${FMP}_allopenchrom.bed
bedtools closest -t first -b loci_tss.bed -a ${FMP}_allopenchromgrp.txt -wa -wb > ${FMP}_allopenchromsnp.tmp
paste <(cut -f 7-12 ${FMP}_allopenchromsnp.tmp) <(cut -f 1-6 ${FMP}_allopenchromsnp.tmp) > ${FMP}_allopenchromsnp.bed


#6. Annotate FINEMAP SNPs for GTEx eQTL in trait-matched tissues
bedtools intersect -a ../../MAURANO-GTEx/${TRAIT}_gtex.bed -b ${FMP}.bed -wa -wb > ${FMP}_gtex.txt
awk 'FNR==NR { split($4,a,","); m[a[1]","a[2]","a[3]]=$4; next } { OFS="\t"; $4=m[$4]; print $0 }' loci_genes.bed ${FMP}_gtex.txt > ${FMP}_gtex.bed

echo -e "snp.chrom\tsnp.chrom.start\tsnp.chrom.end\tsnp.key\tsnp.log10bf\tsnp.strand" > ${FMP}_snp_masterlist.txt
awk 'FNR==NR { m[$10]=$0;next; } $4 in m { print }' <(cat ${FMP}_{gene,tss,func,snpeff,openchrom,allopenchrom,allopenchromsnp,gtex}.bed) ${FMB} | sort | uniq >> ${FMP}_snp_masterlist.txt

# awk 'FNR==NR { m[$4]=$0;next; } $4 in m { print }' <(cat ${FMP}_allopenchromgrp.txt ${FMP}_openchromgrp.txt) ${FMB} > ${FMP}_snp_masterlist.tmp2
# cat ${FMP}_snp_masterlist.tmp1 ${FMP}_snp_masterlist.tmp2 | sort | uniq >> ${FMP}_snp_masterlist.txt

echo -e "snp.key\tsnp.locus\tsnp.name\tgene.key\tgene.gid\tgene.eid\tgene.name\tgene.locus" > ${FMP}_gene_masterlist.txt
awk -F '\t' '$10 != "." && $10 != "" && $4 != "" { OFS="\t"; split($10,b,","); split($4,a,","); split(m[$10],c,"\t"); print $10, b[1], b[2], $4, a[1], a[2], a[3], a[4] }' ${FMP}_{gene,tss,func,snpeff,openchrom,allopenchrom,allopenchromsnp,gtex}.bed | sort | uniq >> ${FMP}_gene_masterlist.txt
