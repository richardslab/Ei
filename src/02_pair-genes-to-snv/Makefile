.DEFAULT_GOAL := all

################################ INPUT PARAMS #################################
# GWAS directory name
trait=glucose
# GWAS sample size
N=375396
# FINEMAP prior stddev
PRIORSTD=0.05
###############################################################################

# note that code below required specific directory layout. If you would like to know more please contact vince.forgetta@ladydavis.ca.

# Default locus size
locussize = 500000
# Merge adjacent loci that are less than value below
mergedist = 50000
# Locus designation (normally supplied as command line parameter).
LOCUS?=L1
# Prefix for output directory.
LPREFIX = regions/$(LOCUS)
# Chromosome containing locus (normally supplied as command line parameter).
CHROM?=1

all: leadSnps.bed leadLoci.bed

# 1. Convert lead SNPs to BED format
leadSnps.bed: ../../../../data/traits/$(trait)/gwas.txt ../../../../data/traits/$(trait)/leadSnps.txt
	awk 'FNR==NR { m[$$1]=$$1; next } ($$1 in m) { print "chr"$$2, $$3, $$3+1, $$1 }' ../../../../data/traits/$(trait)/leadSnps.txt ../../../../data/traits/$(trait)/gwas.txt | cut -b 4- | sort -k 1,1n -k 2,2n | awk '{ print "chr"$$0 }' | uniq -u > leadSnps.bed

# 2. From lead SNPs, output windows of $(locussize) centered on lead SNP to leadLoci.bed, name each locus L#, where # is sequential index sorted by genomic position of the lead SNP.
leadLoci.bed: leadSnps.bed
	sed "s/ /\t/g" leadSnps.bed > leadSnps.tab
	bedtools merge -i leadSnps.tab -d $(mergedist) -o collapse,count -c 4 > leadSnps.merged
    # Pad loci to $(locussize)
	awk -v lsize=$(locussize) '{ c=$$1; e=$$3; s=$$2; r=lsize/2; s=s-r; e=e+r; if(s<0) { s=1 }; printf "%s\t%d\t%d\tL%s\n",c,s,e, NR }' leadSnps.merged > leadLoci.bed
	awk 'c==$$1{ print $$4, $$2-e, $$0, p; e=$$3; p=$$0; c=$$1 }' leadLoci.bed | sort -k2,2nr
#	cut -b 4- leadSnps.bed | sort -k 1,1n -k 2,2n | awk -v lsize=$(locussize) 'BEGIN { OFS="\t" } { fl=(lsize-($$3-$$2))/2; s=$$2-fl; if (s < 0) { s = 0 }; e=$$3+fl; printf "chr%d\t%d\t%d\tL%s\n",$$1,s,e,NR }' > leadLoci.bed


# 3. Prepare input and run FINEMAP
# Generate input for FINEMAP
finemap.1: ../../../../data/traits/$(trait)/gwas.txt leadSnps.bed leadLoci.bed
	mkdir -p regions/
	@echo "Processing "$(LOCUS) $(CHROM) > /dev/stderr
	head -n 1 ../../../../data/traits/$(trait)/gwas.txt > $(LPREFIX).gwas
# Fetch GWAS summary statistics for specified $(LOCUS).
	awk -v LOCUS=$(LOCUS) 'FNR==NR { if ($$4==LOCUS) { sub("chr","",$$1); lc=$$1;ls=$$2;le=$$3; }; next } $$2==lc && $$3>ls && $$3<le { print }' leadLoci.bed ../../../../data/traits/$(trait)/gwas.txt >> $(LPREFIX).gwas

# Fetch genotype data for SNPs at $(LOCUS) based on contents of retreived summary statistics.
finemap.2:
	cut -f 1 -d ' ' $(LPREFIX).gwas | perl -p -e 's/\n/ /g' | perl -p -e 's/\s+$$//' > $(LPREFIX).incl_variants
	bgenix -g $$path/to/UKB/FULL_GENETIC_DATA_RELEASE/v3/imputed/bgen/$(CHROM).bgen -incl-rsids $(LPREFIX).incl_variants > $(LPREFIX).bgen
	bgenix -g $(LPREFIX).bgen -index
	bgenix -g $(LPREFIX).bgen -list > $(LPREFIX).list
	egrep -v '^(#|alternate)' $(LPREFIX).list | cut -f 2 | sort | uniq -u >> $(LPREFIX).uniq
# rm -f $(LPREFIX).bgen $(LPREFIX).bgen.bgi $(LPREFIX).list
	mv $(LPREFIX).bgen $(LPREFIX).orig.bgen 
	mv $(LPREFIX).bgen.bgi $(LPREFIX).orig.bgen.bgi
	mv $(LPREFIX).list $(LPREFIX).orig.list
	awk 'FNR==NR { m[$$1]=$$1; next } $$1 in m { print $$1 }' $(LPREFIX).uniq $(LPREFIX).gwas > $(LPREFIX).extract
	bgenix -g $$path/to/UKB/FULL_GENETIC_DATA_RELEASE/v3/imputed/bgen/$(CHROM).bgen -incl-rsids $(LPREFIX).extract > $(LPREFIX).bgen
	bgenix -g $(LPREFIX).bgen -index
	bgenix -g $(LPREFIX).bgen -list > $(LPREFIX).list
	cut -f 1 $$path/to/UKB/FULL_GENETIC_DATA_RELEASE/v3/imputed/plink.bycroftRand50k/1.fam > $(LPREFIX).incl
	cp $$path/to/UKB/UKBPROJECT/full_release/v3/ukb#####_imp_chr1_v3_s######.sample $(LPREFIX).sample
	ldstore --bgen $(LPREFIX).bgen --bcor $(LPREFIX).bcor --n-threads 8 --samples $$path/to/GENOTYPE_DATA/UKB/BMDandFX-24268/full_release/v3/ukb#####_imp_chr1_v3_s######.sample $(LPREFIX).sample --incl-samples $(LPREFIX).incl
	ldstore --bcor $(LPREFIX).bcor --merge 8
	ldstore --bcor $(LPREFIX).bcor --matrix $(LPREFIX).matrix
	perl -p -e "s/^\s+//g;" $(LPREFIX).matrix | perl -p -e "s/[ ]+/ /g;" > $(LPREFIX).ld

finemap.3:
	echo "rsid chromosome position allele1 allele2 maf beta se" > $(LPREFIX).z
	egrep -v '^(#|alternate)' $(LPREFIX).list | awk '{ print $$2"_"$$6"_"$$7 }' > $(LPREFIX).bgen_alleles
	awk 'FNR==NR { m[$$1]=$$1 ;next } $$1 in m { $$2=sprintf("%02d", $$2); print }' $(LPREFIX).extract $(LPREFIX).gwas | \
	awk 'FNR==NR { m[$$1]=1; next } { a1=$$4; a2=$$5; beta=$$7; maf=$$6; if (maf > 0.5) { maf=1-maf } if ($$1"_"a2"_"a1 in m) { a1=$$5; a2=$$4; beta=-beta }; print $$1,$$2,$$3,a1,a2,maf,beta,$$8 }' $(LPREFIX).bgen_alleles - >> $(LPREFIX).z
# echo "z;bgen;bgi;dose;sample;incl;snp;config;log;n_samples" > $(LPREFIX).in
# echo "$(LOCUS).z;$(LOCUS).bgen;$(LOCUS).bgen.bgi;$(LOCUS).dose;$(LOCUS).sample;$(LOCUS).incl;$(LOCUS).snp;$(LOCUS).config;$(LOCUS).log;$(N)" >> $(LPREFIX).in
	echo "z;ld;bgen;bgi;dose;snp;config;cred;log;n_samples" > $(LPREFIX).in
	echo "$(LOCUS).z;$(LOCUS).ld;;;;$(LOCUS).snp;$(LOCUS).config;$(LOCUS).cred;$(LOCUS).log;$(N)" >> $(LPREFIX).in
	head $(LPREFIX).z

finemap: finemap.1 finemap.2 finemap.3

run_finemap:
	finemap1.3.1 --sss --in-files $(LOCUS).in --n-causal-snps 20 --log --group-snps --corr-config 0.9 --corr-group 0.9 --prior-std $(PRIORSTD) --n-iterations 10000

clean_all:
	rm leadSnps.bed leadLoci.bed regions/*

clean_finemap:
	rm regions/$(LOCUS).*
