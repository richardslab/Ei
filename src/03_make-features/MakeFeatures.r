
#list all annotation files
# all.annot.files <- list.files(paste0(path.main,"annotation_dat/"))
# all.annot.files <- list.files("data/01_annotation-matricies/")
all.annot.files <- list.files("data/01_annotation-matricies/", full.names = TRUE)

#for each trait
for( i.trait in c(1:length(all.annot.files))){
  
dat.file <- all.annot.files[i.trait]
# dat.file <- all.annot.files[1]

raw.annot <- read.table(dat.file,header=TRUE)

trait.name <- unique(raw.annot$trait)
##################################
#creat binary variables for snpeff.impact
raw.annot$snpeff.impact.modifier <- raw.annot$snpeff.impact=="MODIFIER"
raw.annot$snpeff.impact.none <- raw.annot$snpeff.impact=="NONE"
raw.annot$snpeff.impact.low <- raw.annot$snpeff.impact=="LOW"
raw.annot$snpeff.impact.moderate <- raw.annot$snpeff.impact=="MODERATE"
raw.annot$snpeff.impact.high <- raw.annot$snpeff.impact=="HIGH"
##################################
#process GWAS signals
raw.annot$gene.length <- abs(raw.annot$gene.chrom.end - raw.annot$gene.chrom.start) 
raw.annot$snp.tss.dist <- abs(raw.annot$snp.tss.dist)
raw.annot$snp.gene.dist <- abs(raw.annot$snp.gene.dist)
raw.annot$beta <- abs(raw.annot$beta)
raw.annot$z <- abs(raw.annot$z)
##################################

load("src/03_make-features/funcs_collapsing.RData")
snp.collapse <- perpersonsummary(raw.annot, functioncalls.collapses  )

#assign gene names to rows
Genes.names <- table(raw.annot$gene.name)
Genes.names <- Genes.names[Genes.names>0]
rownames(snp.collapse) <- names(Genes.names)
colnames(snp.collapse) <- func.collapses
###########################################################################
###########################################################################
#load true positive genes

omim_genes <- read.table("src/03_make-features/positive-control-omim-genes-hdo-latest.txt",header=TRUE)
drug_genes <- read.table("src/03_make-features/positive-control-drug-targets-knownaction-latest.txt",header=TRUE)
amp_genes <- read.table("src/03_make-features/positive-control-t2d-amp-latest.txt",header=TRUE)

#collect PC genes
trait.genes <- unique( c( as.character( omim_genes$Gene.Symbol[grep(trait.name,omim_genes$Trait)] ) ,
                          as.character( drug_genes$Gene.Symbol[grep(trait.name,drug_genes$Trait)] )  ) )

#For T2D genes we add AMP database
if("t2d" %in% trait.name){
  trait.genes <- unique( c( as.character( omim_genes$Gene.Symbol[grep(trait.name,omim_genes$Trait)] ) ,
                            as.character( drug_genes$Gene.Symbol[grep(trait.name,drug_genes$Trait)] ) ,
                            as.character( amp_genes$Gene.Symbol[grep(trait.name,amp_genes$Trait)] ) 
  ) )
}

#the locus containing PC genes
sel.locus <- unique(raw.annot$locus.name[ raw.annot$gene.name %in% trait.genes ])

#genes in the locus of PC genes
genes.locus <- raw.annot$gene.name[raw.annot$locus.name %in% sel.locus]

#restricted gene data
locus.collapse <- snp.collapse[unique(match(genes.locus, rownames(snp.collapse))),]

#add PC label to each gene
PC.gene <- rep(0,nrow(locus.collapse))
PC.gene[rownames(locus.collapse) %in% trait.genes] <- 1

locus.collapse <- cbind(locus.collapse,PC.gene)

#save the collapsed GWAS genes and PC-locus genes.
save(locus.collapse, snp.collapse , file=paste0("data/02_features/", trait.name,".RData"))

}
