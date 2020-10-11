#new.part.metab <- partition_genes(rownames(comb.metab),  comb.metab, comb.metab[,59] )


partition_genes <- function(gene.names,  X, Y.response ){
  
#unique gene IDs

unique.names <- unique(gene.names)

set.seed(99)

uniq.gene.entry <- NULL
uniq.y <- NULL
uniq.id <- NULL
for(temp.gene.id in unique.names){
  
  #print(temp.gene.id)
  
  GWAS.per.gene <-  which( gene.names  == temp.gene.id)
  
  if( length(GWAS.per.gene)==1 ){
    rand.gene.id <- GWAS.per.gene
    temp.gene.entry <- X[GWAS.per.gene,]
    temp.uniq.y <- Y.response[GWAS.per.gene]
  }
  
  #if multiple rows
  else{
    
    #if all negative genes
    if( sum(Y.response[GWAS.per.gene])==0 ){
    rand.gene.id <- sample(GWAS.per.gene,1 )
    temp.gene.entry <- X[rand.gene.id,]
    temp.uniq.y <- Y.response[rand.gene.id]
    }
    
    #if at least 1 positive gene
    else{
      
      #if only 1 positive
      if( length(GWAS.per.gene[Y.response[GWAS.per.gene]==1])==1 )rand.gene.id <- GWAS.per.gene[Y.response[GWAS.per.gene]==1]
      #if multiple positives
      else{   rand.gene.id <- sample( GWAS.per.gene[Y.response[GWAS.per.gene]==1] ,1 ) }
      
      temp.gene.entry <- X[rand.gene.id,]
      temp.uniq.y <- Y.response[rand.gene.id]
      
    }
    
    
    
  }
  
  uniq.gene.entry <- rbind(uniq.gene.entry, temp.gene.entry)
  uniq.y <- c(uniq.y, temp.uniq.y)
  uniq.id <- c( uniq.id, rand.gene.id)
}


  
  list(uniq.names = unique.names, uniq.x = uniq.gene.entry, uniq.y = uniq.y, uniq.id = uniq.id )
  
  

  
  
}









######################################################
######################################################
######################################################
######################################################
######################################################
#testing the function

if(FALSE){
  comb.metab
  
  #the possibly repeated gene names  
  gene.names <- rownames(comb.metab)
  
  chars.name = "PC_exp"
  f <-  as.data.frame(comb.metab)
  X <-  as.data.frame(f)
  
  X <- comb.metab
  #the response value 0 or 1
  id.response <- grep(chars.name,colnames(comb.metab))
  Y.response <- X[,id.response]
  
  res.part <- partition_genes(gene.names,  X, Y.response )
  
}

