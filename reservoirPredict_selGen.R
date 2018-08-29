# Babayan, Orton & Streicker 
# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes 
# Reservoir host prediction from selected genomic features
# https://www.h2o.ai/products/h2o/ 

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o)
library(dplyr) 
library(reshape2)
library(matrixStats)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Start h2o JVM
localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# Read data from file
f1<-read.csv(file="BabayanEtAl_VirusData.csv",header=T)
fis<-read.csv(file="featureImportance_reservoir.csv",header=T)

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)

# Feature selection (simplify dataset to required columns)
nfeats<-50
totalfeats<-length(fis$vimean)
f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
gen.feats<-as.character(fis$X[f])
f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector",gen.feats)] 

# Remove orphans
f2<-subset(f1,f1$Reservoir!="Orphan")
f<-droplevels(f2)
o<-subset(f1,f1$Reservoir=="Orphan")
orphans<-o[with(o,order(Viral.group,Virus.name,decreasing=T)),]

# Group selection based on thresholds
t<-15 # threshold for minimum sample size of groups
s<-.7 # proportion in the training set
host.counts<-table(f$Reservoir)
min.t<-host.counts[host.counts>=t] # minimum number of viruses per host group
f_st3<-f[f$Reservoir %in% c(names(min.t)),]
f_st3<-droplevels(f_st3)
f_st3$SeqName2<-do.call(rbind,strsplit(as.character(f_st3$Genbank.accession),"[.]"))[,1]

# Rare hosts
rare<-f[!f$Reservoir %in% c(names(min.t)),]
rare<-droplevels(rare)
rare$SeqName2<-do.call(rbind,strsplit(as.character(rare$Genbank.accession),"[.]"))[,1]

# Number and names of host taxa
ntax<-length(unique(f_st3$Reservoir))
bp<-as.character(sort(unique(f_st3$Reservoir)))

# Sample split of training/test to get counts in each 
trains<-f_st3 %>% group_by(Reservoir) %>%
  filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
testval<-subset(f_st3,!(f_st3$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set

optims<-testval %>% group_by(Reservoir) %>%
  filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
tests<-subset(testval,!(testval$Genbank.accession %in% optims$Genbank.accession)) # ref numbers in testval set absent from test set    
ntest<-dim(tests)[1]

# Remove unneeded files
rm(f,f1,f2,fis)

# Train many models
set.seed(78910) 
nloops<-550
lr<-c()
md<-c()
sr<-c()
csr<-c()
nt<-c()
mr<-c()
accuracy.st3<-c()

pc.accuracy<-matrix(nrow=nloops,ncol=ntax)
test.record<-matrix(nrow=ntest,ncol=nloops)
nfeatures<-length(gen.feats)
vimps<-matrix(nrow=nfeatures,ncol=nloops)

for (i in 1:nloops){
    # Stratified random
    trains<-f_st3 %>% group_by(Reservoir) %>%
      filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
    testval<-subset(f_st3,!(f_st3$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set
  
    optims<-testval %>% group_by(Reservoir) %>%
      filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
    tests<-subset(testval,!(testval$Genbank.accession %in% optims$Genbank.accession)) # ref numbers in testval set absent from test set    
  
    trains<-droplevels(trains)
    tests<-droplevels(tests)
    optims<-droplevels(optims)
    test.record[,i]<-as.character(tests$Genbank.accession)
    
    set<-c("Reservoir",gen.feats) 
    
    # Build training, optimization, validation, rare and orphan datasets
    f1_train<-trains[,c(set)] # this is the full training dataset with genomic features and known reservoir associations

    testID<-tests$Virus.name
    f1_test<-tests[,c(set)] 

    optID<-optims$Virus.name
    f1_opt<-optims[,c(set)] 
    
    set<-c(gen.feats) # only genomic features for orphan and rare viruses
    f1_orphan<-orphans[,c(set)] 
    f1_rare<-rare[,c(set)] 
    
    # Convert to h2o data frames    
    train<-as.h2o(f1_train)
    test<-as.h2o(f1_test)
    opt<-as.h2o(f1_opt)
    orp<-as.h2o(f1_orphan)
    rar<-as.h2o(f1_rare)
    
    # Cleanup
    rm(f1_train,f1_test,f1_opt,f1_orphan,f1_rare)
    
    # Identity the response column
    y <- "Reservoir"

    # Identify the predictor columns
    x <- setdiff(names(train), y)

    # Convert response to factor
    train[,y] <- as.factor(train[,y])
    test[,y] <- as.factor(test[,y])
    opt[,y] <- as.factor(opt[,y])

    # Train and validate a grid of GBMs
    gbm_params <- list(learn_rate = c(.001,seq(0.01, 0.2, .02)),
    max_depth = seq(6, 15, 1),
    sample_rate = seq(0.6, 1.0, 0.1),
    col_sample_rate = seq(0.5, 1.0, 0.1),
    ntrees=c(100,150,200),
    min_rows=c(5,8,10))

    search_criteria <- list(strategy = "RandomDiscrete",
    max_models = 500,
    stopping_rounds=10,
    stopping_metric="misclassification",
    stopping_tolerance=1e-3)

    gbm_grid <- h2o.grid("gbm", x = x, y = y,
    grid_id = "gbm_grid",
    training_frame = train,
    validation_frame = opt, 
    seed = 1,
    hyper_params = gbm_params,
    search_criteria = search_criteria)

    gbm_gridperf <- h2o.getGrid(grid_id = "gbm_grid",
    sort_by = "accuracy",
    decreasing = TRUE)

    # Grab the model_id for the top GBM model
    best_gbm_model_id <- gbm_gridperf@model_ids[[1]]
    best_gbm <- h2o.getModel(best_gbm_model_id)
    perf <- h2o.performance(best_gbm, test)
    
    # Record best settings
    lr[i]<-as.numeric(gbm_gridperf@summary_table[1,2]) # learn_rate
    sr[i]<-as.numeric(gbm_gridperf@summary_table[1,6]) # sample_rate
    md[i]<-as.numeric(gbm_gridperf@summary_table[1,3]) # maxdepth
    csr[i]<-as.numeric(gbm_gridperf@summary_table[1,1]) # col_sample_rate
    mr[i]<-as.numeric(gbm_gridperf@summary_table[1,4]) # col_sample_rate
    nt[i]<-as.numeric(gbm_gridperf@summary_table[1,5]) # col_sample_rate
    
    # Extract confusion matrix and calculate accuracies
    cm1<-h2o.confusionMatrix(perf)
    nclass<-length(unique(trains$Reservoir))
    cm2<-cm1[1:nclass,1:nclass]
    cm<-as.matrix(cm2)

    norm_cm<-cm/rowSums(cm)
    accuracy.st3[i]=sum(diag(cm))/sum(cm)
    pc.accuracy[i,]<-t(diag(cm)/rowSums(cm))
    
    write.csv(norm_cm,file=paste("Reservoir_selGen_CM",i,".csv"))
    
    # Retreive feature importance
    vi <- h2o.varimp(best_gbm)
    data2  <- vi[order(vi[,1],decreasing=FALSE),]
    vimps[,i]<-data2[,4]

    # Orphan predictions
    orp.pred <- h2o.predict(best_gbm, orp)
    df<-orp.pred[,c(2:(ntax+1))]
    df2<-as.data.frame(df)
    row.names(df2)<-orphans$Virus.name
    write.csv(df2,file=paste("Orphans",i,".csv",sep="_")) 
    
    # Rare predictions
    rar.pred <- h2o.predict(best_gbm, rar)
    df<-rar.pred[,c(2:(ntax+1))]
    df2<-as.data.frame(df)
    row.names(df2)<-rare$Virus.name
    write.csv(df2,file=paste("ST5_RareVirus",i,".csv",sep="_")) 
    
    # Test set predictions
    test.pred<-h2o.predict(best_gbm,test[,2:length(names(test))]) 
    df2<-as.data.frame(test.pred)
    row.names(df2)<-testID
    write.csv(df2,file=paste("TestPred",i,".csv",sep="_"))  
    
    # Clean up
    h2o.rm("gbm_grid")
    rm(gbm_grid,best_gbm,train,test,opt,df2,optims)
}

accs<-data.frame(accuracy.st3,pc.accuracy,lr,sr,md,csr)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(vimps,file="Reservoir_SelGen50_FI.csv",row.names = T)
write.csv(accs,file="Reservoir_SelGen50_out.csv",row.names=F)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
vecs<-table(tests$response)
chanceAccurate<-round(sum(prob*vecs),digits=0)
tot<-sum(vecs)
nullAcc<-chanceAccurate/tot
print(nullAcc)