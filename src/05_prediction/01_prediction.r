#we already predicted each of the 12 traits using data from the remaining traits
#in the xgboost training and test phase.
#here we just collect and aggregate the predicted probabilities

all.locus.prob <- NULL
all.locus.y <- NULL
names.genes <- NULL
all.trait <- NULL


#the folder containing the trained xgboost model
path.load <- "~"

for(i.trait in c(1:12)){

  load(file=paste0(path.load,i.trait,".RData"),verbose=TRUE)
  
  #the name of tested genes
  names.genes <- c(names.genes, rownames(test.x)[rownames(test.x)%in%test.locus.genes])
  #predicted probability to be PC
  locus.pred.y <- fit[rownames(test.x)%in%test.locus.genes]
  #true label of PC
  locus.true.y <- test.y[rownames(test.x)%in%test.locus.genes]
  
  all.locus.prob <- c(all.locus.prob, locus.pred.y)
  all.locus.y    <- c(all.locus.y,    locus.true.y)
  all.trait      <- c(all.trait, rep(i.trait,length(locus.true.y)))
  
  #PR curve
  pred.prob <- locus.pred.y
  fg <- pred.prob[locus.true.y == 1]
  bg <- pred.prob[locus.true.y == 0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T,rand.compute=T)
  
  #ROC curve
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  #if plot ROC/PR curve
  if(TRUE){
    
    trait.name <- substr(traits.files, 1, 3)[i.trait]
    if(trait.name=="tg.")trait.name <- "tg_"
    
    setwd(path.save)
    
    jpeg(paste0(trait.name,"-PRcurve.jpeg"))
    plot(pr,main=paste0(trait.name,"-","PR"))
    dev.off()
    
    jpeg(paste0(trait.name,"-ROCcurve.jpeg"))
    plot(roc,main=paste0(trait.name,"-","PR"))
    dev.off()
    
    
  }
  
}

pred.tab <- data.frame(
  all.locus.prob ,
  all.locus.y ,
  names.genes ,
  traits= gsub( ".RData", "",traits.files[all.trait] ) 
)

