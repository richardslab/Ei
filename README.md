# Effector Index Algorithm

The algorithm consists of 5 steps.

1. Statistical fine-mapping (`src/01_finemapping`): GWAS summary statistics are hamonized (`01_format-gwas-sumstats.sh`) then LD clumped to determine lead SNPs (`02_LD-clump-sumstats.sh`). Loci are determined from the lead SNPs and FINEMAP is run for GWAS summary statistics for each of the loci (03_run-finemap.sh).

2. Annotating causal SNVs and linking to nearby genes (`src/02_pair-genes-to-snv`): FINEMAP results are filtered for loci that failed convergence or overlap the MHC. FINEMAP SNVs with log10(BF) > 2 were retained (`01_clean-finemap.sh`). Causal SNVs are annotated for various functional annotations and genes at the locus are paired to one or more SNVs(`02_annotate.sh`), resulting in an annotation matrix (`03_build-annotation-matrix.sh`).

3. Building features from the annotation matrix (`src/03_make-features`): The annotation matrix  may contain multiple SNVs paired to a given gene and thus multiple annotation values (such as coding impact, GWAS effect size, etc) may be  associated to each gene. To resolve each annotation to a single value per gene (hereby referred to as features), multiple summarizations were conducted (sum, min, mean, max) as well as scaling of these values based on their distance to the paired gene (`01_make-features.r`).

4. Training prediction models (`src/04_train-models`): First, perform unique partition on genes across 12 traits before model fitting to avoid using the same gene multiple times for different traits thus remove overfitting bias (01_partition.r). Model training was conducted using XGBoost in a leave-one-out approach, whereby each trait-specific model was trained using data from the remaining traits (`02_trainModel.r`).
5. Prediction of the Effector Index for a GWAS trait (`src/05_prediction`):  Extract predictions from the left-out trait (01_prediction.r).


