library(data.table)
rm(list = ls(all.names = TRUE))

args <- commandArgs(trailingOnly = TRUE)

## GWAS trait
trait.dir <- "M01494"
#trait.name <- "M01494"
#trait.dir <- args[1]
#trait.name <- args[2]

## File prefix, usually denotes the minimum FINEMAP log10bf.
FMP <- "fm2" 

get.gene.pos <- function(gene.name, gene.chrom.start, gene.chrom.end, gene.strand, snp.chrom.start) {
## Calculate position of SNP relative to gene boundaries. 
## * Upstream distance is reported from the TSS and is negative.
## * Downstream distance is reported from the TES and is positive.
## * Gene body overlap is reported as 0.
  gene.pos <- NA
  if (gene.strand == "+"){
    if (snp.chrom.start < gene.chrom.start){
      gene.pos <- snp.chrom.start - gene.chrom.start
    }else{
      gene.pos <- snp.chrom.start - gene.chrom.end
    }
  }
  
  if (gene.strand == "-"){
    if (snp.chrom.start < gene.chrom.start){
      gene.pos <- gene.chrom.start - snp.chrom.start
    }else{
      gene.pos <- gene.chrom.end - snp.chrom.start
    }
  }
  if (snp.chrom.start >= gene.chrom.start && snp.chrom.start <= gene.chrom.end){
    gene.pos <- 0L
  }
  return(gene.pos)
}

## Load unique list of SNPs present across all annotations
DT_snp_master <- fread(paste0(trait.dir, "/annotation/",FMP, "_snp_masterlist.txt"), header = TRUE)

## Load unique list of genes present across all annotations
DT_gene_master <- fread(paste0(trait.dir, "/annotation/",FMP, "_gene_masterlist.txt"), header = TRUE)
## Merge the two master lists and generate unique key for each SNP-gene pair
DT_master <- merge(DT_gene_master, DT_snp_master, by="snp.key", all.x=TRUE)
DT_master[, key := paste0(snp.key,":",gene.key)]
message("0. Merge keys for SNPs and genes, row total = ", DT_master[,.N])

DT_leadsnps <- fread(paste0(trait.dir, "/annotation/leadSnps.bed"), header=FALSE, col.names = c("chrom", "chrom.start","chrom.end", "snp.name"))
DT_master[, is.leadsnp := snp.name %in% DT_leadsnps$snp.name]

# rm("DT_snp_master", "DT_gene_master")

# Gene information
DT_gene <- fread(paste0(trait.dir, "/annotation/","loci_gene_overlaps.bed"),
              col.names = c("locus.chrom", "locus.chrom.start", "locus.chrom.end", "locus.name",
                            "gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                            "gene.score", "gene.strand"))
DT_gene[, gene.key := paste0(gene.key,",",locus.name)]
DT_master <- merge(DT_master, DT_gene, by="gene.key", all.x=TRUE)
message("1. Merge gene information, row total = ", DT_master[,.N])
DT_master[, gene.tss := ifelse(gene.strand == "+", gene.chrom.start, gene.chrom.end)]
DT_master[, snp.tss.dist := ifelse(gene.strand == "+", snp.chrom.start- gene.tss, gene.tss - snp.chrom.start)]

DT_master[, snp.gene.dist := get.gene.pos(gene.name, gene.chrom.start, gene.chrom.end, gene.strand, snp.chrom.start), by=seq_len(nrow(DT_master))]

DT_fm_gene <- fread(paste0(trait.dir, "/annotation/",FMP, "_gene.bed"), 
                    col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                                  "gene.score", "gene.strand",
                                  "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                  "snp.log10bf", "snp.strand"
                                  ))

DT_fm_gene[, key := paste0(snp.key,":",gene.key)]
DT_master[, is.gene.overlap := key %in% DT_fm_gene$key]
message("2. Annotating SNP/gene overlaps, row total = ", DT_master[,.N])

DT_fm_tss <- fread(paste0(trait.dir, "/annotation/",FMP, "_tss.bed"),
                col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                              "gene.score", "gene.strand",
                              "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                              "snp.log10bf", "snp.strand", "snp.tss.dist"
                              ))
DT_fm_tss[, key := paste0(snp.key,":",gene.key)]
DT_master[, is.nearest.to.tss := key %in% DT_fm_tss$key]
# DT_master <- merge(DT_master, DT_fm_tss[, .(key, snp.tss.dist)], by="key", all.x=TRUE)
message("3. Annotating SNP nearest to each gene TSS, row total = ", DT_master[,.N])

DT_finemap <- fread(paste0(trait.dir, "/annotation/",FMP, ".txt"), header = TRUE)
DT_finemap[, snp.key := paste0(locus,",",rsid)]
DT_master <- merge(DT_master, DT_finemap[, .(snp.key, maf, beta, se, z, prob, log10bf, log10bf_group, prob_group)], by="snp.key", all.x=TRUE)


DT_fm_func <- fread(paste0(trait.dir, "/annotation/",FMP, "_func.bed"),
                   col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                                 "gene.score", "gene.strand",
                                 "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                 "dbsnp.func", "snp.strand"
                   ))
DT_fm_func[, key := paste0(snp.key,":",gene.key)]
DT_master <- merge(DT_master, DT_fm_func[,.(key, dbsnp.func)], by="key", all.x=TRUE)
message("4. Annotating dbSNP functional consequence, row total = ", DT_master[,.N])
DT_master[is.na(dbsnp.func), dbsnp.func := "none"]
DT_master[, is.dbsnp.delit := dbsnp.func %in% c("missense", "stop-loss","nonsense","frameshift", "cds-indel", "splice-5", "splice-3")]
DT_master[is.na(is.dbsnp.delit), is.dbsnp.delit := FALSE]

DT_fm_snpeff <- fread(paste0(trait.dir, "/annotation/",FMP, "_snpeff_distinct.bed"),
                    col.names = c("chrom", "chrom.start", "chrom.end", "key", 
                                  "snpeff.impact"))

DT_master <- merge(DT_master, DT_fm_snpeff[,.(key, snpeff.impact)], by="key", all.x=TRUE)
message("5. Annotating SNPEff functional impact, row total = ", DT_master[,.N])
DT_master[is.na(snpeff.impact), snpeff.impact := "NONE"]

# Code has bug, returns NA for SNPs with multiple impacts. Code below fixes this, returning integer for highest impact.
# DT_master[, snpeff.rank := as.numeric(factor(snpeff.impact, 
#                                       levels=c("HIGH", "MODERATE", "LOW", "MODIFIER", "NONE"),
#                                       ordered=TRUE))]


DT_master[ snpeff.impact %like% "NONE", snpeff.rank := 0]
DT_master[ snpeff.impact %like% "MODIFIER", snpeff.rank := 1]
DT_master[ snpeff.impact %like% "LOW", snpeff.rank := 2]
DT_master[ snpeff.impact %like% "MODERATE", snpeff.rank := 3]
DT_master[ snpeff.impact %like% "HIGH", snpeff.rank := 4]

DT_master[, is.snpeff.delit := snpeff.impact %like% "(HIGH|MODERATE)"]
DT_master[is.na(is.snpeff.delit), is.snpeff.delit := FALSE]

#DT_fm_openchomgrp <- fread(paste0(trait.dir, "/annotation/",FMP, "_openchromgrp.bed"),
#                      col.names = c("snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
#                                    "openchrom.count", "snp.strand"
#                      ))
#DT_master <- merge(DT_master, DT_fm_openchomgrp[,.(snp.key, openchrom.count)], by="snp.key", all.x=TRUE)

# CLOSEST TRAIT-MATCHED OPENCHROM SNP TO EACH GENE
DT_fm_openchom <- fread(paste0(trait.dir, "/annotation/",FMP, "_openchromgrp.txt"),
                        col.names = c("snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                      "snp.in.trait.DHS", "snp.strand"
                        ))
if (nrow(DT_fm_openchom) > 0){
  DT_master <- merge(DT_master, DT_fm_openchom[, .(snp.key, snp.in.trait.DHS)], by="snp.key", all.x=TRUE)
  DT_master[is.na(snp.in.trait.DHS), snp.in.trait.DHS := 0]
  DT_tmp <- DT_master[snp.in.trait.DHS > 0, min(abs(snp.tss.dist)), by=.(gene.name)]
  DT_tmp[, key := paste0(gene.name, ",", V1)]
  DT_master[, nearest.trait.DHS.from.gene := paste0(gene.name, ",", abs(snp.tss.dist)) %in% DT_tmp$key]
} else{
  DT_master[, snp.in.trait.DHS := 0]
  DT_master[, nearest.trait.DHS.from.gene := FALSE]
}
# CLOSEST GENE FROM A TRAIT-MATCHED OPENCHROM SNP
DT_fm_openchom <- fread(paste0(trait.dir, "/annotation/",FMP, "_openchromsnp.bed"),
                        col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                                      "gene.score", "gene.strand",
                                      "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                      "nearest.gene.from.trait.DHS", "snp.strand"
                        ))
if (nrow(DT_fm_openchom) > 0){
  DT_fm_openchom[, key := paste0(snp.key,":",gene.key)]
  DT_master[, nearest.gene.from.trait.DHS := key %in% DT_fm_openchom$key]
} else {
  DT_master[, nearest.gene.from.trait.DHS := FALSE]
}
# CLOSEST OPENCHROM SNP TO EACH GENE
DT_fm_openchom_all <- fread(paste0(trait.dir, "/annotation/",FMP, "_allopenchromgrp.txt"),
                        col.names = c("snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                      "snp.in.DHS", "snp.strand"
                        ))
DT_master <- merge(DT_master, DT_fm_openchom_all[, .(snp.key, snp.in.DHS)], by="snp.key", all.x=TRUE)
DT_master[is.na(snp.in.DHS), snp.in.DHS := 0]
DT_tmp <- DT_master[snp.in.DHS > 0, min(abs(snp.tss.dist)), by=.(gene.name)]
DT_tmp[, key := paste0(gene.name, ",", V1)]
DT_master[, nearest.DHS.from.gene := paste0(gene.name, ",", abs(snp.tss.dist)) %in% DT_tmp$key]

# CLOSEST GENE FROM A TRAIT-MATCHED OPENCHROM SNP
message("6. Annotating openchrom, row total = ", DT_master[,.N])
DT_fm_openchom <- fread(paste0(trait.dir, "/annotation/",FMP, "_allopenchromsnp.bed"),
                        col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                                      "gene.score", "gene.strand",
                                      "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                      "nearest.gene.from.DHS", "snp.strand"
                        ))
DT_fm_openchom[, key := paste0(snp.key,":",gene.key)]
DT_master[, nearest.gene.from.DHS := key %in% DT_fm_openchom$key]



#DT_tmp <- DT_master[snp.in.DHS > 0, .(min(abs(snp.tss.dist))), by=.(snp.key)]
#DT_tmp[, key := paste0(snp.key, ",", V1)]
#DT_master[, nearest.gene.from.DHS := paste0(snp.key, ",", abs(snp.tss.dist)) %in% DT_tmp$key]

#DT_fm_openchom <- fread(paste0(trait.dir, "/annotation/",FMP, "_openchrom.bed"),
#                        col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
#                                      "gene.score", "gene.strand",
#                                      "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
#                                      "trait.openchrom.count", "snp.strand"
#                        ))
#DT_fm_openchom[, key := paste0(snp.key,":",gene.key)]
#DT_master <- merge(DT_master, DT_fm_openchom[!duplicated(snp.key),.(snp.key, trait.openchrom.count)], by="snp.key", all.x=TRUE)
#message("6. Annotating trait openchromatin peak with SNP that is nearest to each gene TSS, row total = ", DT_master[,.N])
#DT_master[is.na(trait.openchrom.count), trait.openchrom.count := 0]
#DT_master[, is.nearest.trait.openchrom := key %in% DT_fm_openchom$key]

#DT_fm_openchom_all <- fread(paste0(trait.dir, "/annotation/",FMP, "_allopenchrom.bed"),
#                        col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
#                                      "gene.score", "gene.strand",
#                                      "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
#                                      "all.openchrom.count", "snp.strand"
#                        ))
#DT_fm_openchom_all[, key := paste0(snp.key,":",gene.key)]
#DT_master <- merge(DT_master, DT_fm_openchom_all[!duplicated(snp.key),.(snp.key, all.openchrom.count)], by="snp.key", all.x=TRUE)
#message("7. Annotating all openchromatin peak with SNP that is nearest to each gene TSS, row total = ", DT_master[,.N])
#DT_master[is.na(all.openchrom.count), all.openchrom.count := 0]
#DT_master[, is.nearest.all.openchrom := key %in% DT_fm_openchom$key]

# GTEx eQTL for SNP-GENE pairs
message("7. Annotating GTEx, row total = ", DT_master[,.N])
DT_fm_gtex <- fread(paste0(trait.dir, "/annotation/",FMP, "_gtex.bed"),
                        col.names = c("gene.chrom", "gene.chrom.start", "gene.chrom.end", "gene.key", 
                                      "gene.score", "gene.strand",
                                      "snp.chrom", "snp.chrom.start", "snp.chrom.end", "snp.key", 
                                      "gtex.score", "snp.strand"
                        ))
if (nrow(DT_fm_gtex) > 0){
  DT_fm_gtex[, key := paste0(snp.key,":",gene.key)]
  DT_master[, in.gtex := key %in% DT_fm_gtex$key]
}else{
  DT_master[, in.gtex := FALSE]
}

#DT_fm_gtex$gene.gid <- unlist(strsplit(DT_fm_gtex$gene.ids, split = ","))[seq(1,nrow(DT_fm_gtex)*3,by=3)]
#DT_fm_gtex[, gtex.key := paste0(snp.key,":",gene.gid)]
#DT_master[, gtex.key := paste0(snp.key,":",gene.gid)]
#DT_master[, in.gtex := gtex.key %in% DT_fm_gtex$key]
message("Duplicated keys (must be 0) = ", DT_master[duplicated(key), .N])

if (TRUE) {
  DT_master[, key := NULL]
  DT_master[, gene.locus := NULL]
  DT_master[, snp.chrom := NULL]
  DT_master[, snp.chrom.end := NULL]
  DT_master[, snp.strand := NULL]
  DT_master[, snp.log10bf := NULL]
  DT_master[, gene.chrom := NULL]
  DT_master[, gene.score := NULL]
  DT_master[, snp.key := NULL]
  DT_master[, gene.key := NULL]
  DT_master[, snp.pos := snp.chrom.start]
  DT_master[, snp.chrom.start := NULL]
}

DT_master[, trait := trait.name]

setcolorder(DT_master, c('trait',
'snp.locus',
'snp.name',
'gene.gid',
'gene.eid',
'gene.name',
'snp.pos',
'is.leadsnp',
'locus.chrom',
'locus.chrom.start',
'locus.chrom.end',
'locus.name',
'gene.chrom.start',
'gene.chrom.end',
'gene.strand',
'gene.tss',
'snp.tss.dist',
'snp.gene.dist',
'is.gene.overlap',
'is.nearest.to.tss',
'maf',
'beta',
'se',
'z',
'prob',
'log10bf',
'log10bf_group',
'prob_group',
'dbsnp.func',
'is.dbsnp.delit',
'snpeff.impact',
'snpeff.rank',
'is.snpeff.delit',
'snp.in.trait.DHS',
'nearest.trait.DHS.from.gene',
'nearest.gene.from.trait.DHS',
'snp.in.DHS',
'nearest.DHS.from.gene',
'nearest.gene.from.DHS',
'in.gtex'))

write.table(DT_master, file=paste0(trait.dir,"/annotation/",trait.dir,"-annotation-matrix.txt"), quote=FALSE, sep="\t", row.names = FALSE)




