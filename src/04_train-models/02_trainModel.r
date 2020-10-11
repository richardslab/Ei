
#For a single trait i.trait = 1, 2,...,12
i.trait <- 2

#load the partitioned data

load("data/03_partitioned-gene-sets/GenPartitionLocus.RData",verbose=TRUE) 

#define training traits and test traits
id.tr <- which(new.new.all.genes$GWAS!=i.trait)
id.ts <- which(new.new.all.genes$GWAS==i.trait)
test.locus.genes <- new.new.all.genes$gene.names[id.ts]

#QC of training datasets
#remove gene names and sd=0 features.
unique.genes <- new.new.all.genes[,-grep("gene.names",colnames(new.new.all.genes))]
unique.genes<- unique.genes[,apply(unique.genes,2,sd)!=0]

#now test if there are remaining highly correlated features (>0.99)
tmp <- cor(unique.genes);tmp[upper.tri(tmp)] <- 0;diag(tmp) <- 0
if(sum(apply(tmp,2,function(x) any(abs(x) > 0.99)) )>0)unique.genes <- unique.genes[,-which( apply(tmp,2,function(x) any(abs(x) > 0.99)) == 1 )]

##########################################################


# path.work <- "/scratch/greenwood/lai.jiang/BACKUP/lai.jiang/Brent/Github/Submission/save_process/"
path.work <- "data/02_features/"
traits.files <- list.files(path.work)

load( paste0( path.work, traits.files[i.trait]),verbose=TRUE)

##########################################################
#Note that all the summarized features start with fn.as.
train.y <- unique.genes[id.tr,grep("PC.gene",colnames(unique.genes))]
train.y <- as.numeric(train.y>0)

test.y <- unique.genes[id.ts,grep("PC.gene",colnames(unique.genes))]
test.y <- as.numeric(test.y>0)

#remove the variable PC.gene and GWAS id from regressors
fsel <- grep("fn.",colnames(unique.genes))

train.x <- unique.genes[id.tr, fsel ]
test.x  <- unique.genes[id.ts, fsel ]

#standardization
mn.x <- apply(train.x, 2,mean)
sd.x <- apply(train.x,2, sd)

#if sd.x==0, we do not need to perform standardization for this feature
if(sum(sd.x==0)>0){
  mn.x[sd.x==0] <- 0
  sd.x[sd.x==0] <- 1
}
#standardization
sd.tr.x <- sweep(train.x, 2, mn.x, FUN = "-")
sd.tr.x <- sweep(sd.tr.x, 2, sd.x, FUN = "/")
sd.ts.x <- sweep(test.x, 2, mn.x, FUN = "-")
sd.ts.x <- sweep(sd.ts.x, 2, sd.x, FUN = "/")

train.x <- sd.tr.x
test.x  <- sd.ts.x
########################################################################
########################################################################
########################################################################
#train xgboost model
library(xgboost)
library(mlr)
library(pROC)

trainset<-data.frame(cbind(train.x,train.y))
trainset$train.y<-as.factor(trainset$train.y)
colnames(trainset)<-c(colnames(train.x),"Response")


lrn <- makeLearner("classif.xgboost",predict.type = "prob")
lrn$par.vals <- list(objective="binary:logistic", eval_metric="aucpr", nrounds=150L, eta=0.1)
params <- makeParamSet(makeDiscreteParam("booster",values = "gbtree"), 
                       makeIntegerParam("max_depth",lower = 3L,upper = 10L),
                       makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                       makeNumericParam("colsample_bytree",lower = 0.5,upper = 1),
                       makeNumericParam("lambda",lower = 0,upper = 5),
                       makeNumericParam("gamma",lower = 1,upper = 10)
)

rdesc <- makeResampleDesc("CV",stratify = T,iters=3L)
ctrl <- makeTuneControlRandom(maxit =10L)
library(parallel)
library(parallelMap) 
#parallelStartSocket(cpus =detectCores())
######################################################################
#we assign weights to prioritize the identification of positive cases
weight <- rep(1,length(train.y))
weight[which(train.y==0)] <- (sum(train.y==1))/(sum(train.y==0))

#training phase of xgboost
traintask <- makeClassifTask (data=trainset,target = "Response",weights=weight)
mytune <- tuneParams(learner = lrn, task = traintask, resampling = rdesc, measures = mlr::gpr, par.set = params, control = ctrl, show.info = T)
lrn_tune <- setHyperPars(lrn,par.vals = mytune$x)
xgmodel <- mlr::train(learner = lrn_tune,task = traintask)


#test set specification and prediction
testset<-data.frame(cbind(test.x,test.y))
testset$test.y<-as.factor(testset$test.y)
colnames(testset)<-c(colnames(test.x),"Response")
testtask <- makeClassifTask (data = testset,target = "Response")
xgpred <- predict(xgmodel,testtask)
fit <- xgpred[["data"]][["prob.1"]]

#performance
roc_obj <- pROC::roc( test.y,  fit )
rec.auc <- pROC::auc(roc_obj)
#feature weights
FI.res<-mlr::getFeatureImportance(xgmodel)
fts.weight <- as.numeric(FI.res$res)

save( xgmodel, rec.auc, fit, id.tr, id.ts, 
     lrn_tune, traintask,
     train.x, test.x,      mn.x,     sd.x, 
     train.y, test.y,
     test.locus.genes,
     file=paste0("data/04_models/", sprintf("%d.RData",i.trait)))

