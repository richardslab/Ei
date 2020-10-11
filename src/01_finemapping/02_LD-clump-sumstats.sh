## Reference panel consists of random subset of 50k individuals from the White-British subset as defined by Bycroft et al. https://www.nature.com/articles/s41586-018-0579-z.
## Imputed genotype data was converted to PLINK format using the method below
## cut -f 2 path/to/UKB/FULL_GENETIC_DATA_RELEASE/v3/imputed/bgen.stats/ukb_mfi_chr${PBS_ARRAYID}_v3.txt | sort | uniq -d > ${PBS_ARRAYID}.dups
## plink2 --bgen path/to/UKB/FULL_GENETIC_DATA_RELEASE/v3/imputed/bgen/${PBS_ARRAYID}.bgen \
##    --sample path/to/UKB/UKBPROJECT/full_release/v3/ukb#####_imp_chr1_v3_s######.sample \
##    --keep rand.50k.FID_IID.txt \
##    --make-bed \
##    --out ${PBS_ARRAYID} \
##    --threads 8 \
##    --memory 40000 \
##    --exclude ${PBS_ARRAYID}.dups

REFP_DIR=path/to/plink.bycroftRand50k

for i in {1..22};
do
    plink --clump gwas.assoc \
        --clump-p1 5e-8 \
        --clump-p2 5e-4 \
        --clump-r2 0.01 \
        --clump-kb 250 \
        --bfile ${REFP_DIR}/${i} --out ${i}
done

awk '$0 ~ /\S+/ && $0 !~ /SNP/ { print $3 }' *.clumped > leadSnps.txt


