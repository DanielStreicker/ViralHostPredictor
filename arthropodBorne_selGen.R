"""
Babayan, Orton & Streicker

Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

-- Arthropod-bone transmission prediction from selected genomic features
"""

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o
library(dplyr)
library(reshape2)
library(matrixStats)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Start h2o JVM
localh20<-h2o.init(nthreads = -1)

# Read data from file
f1<-read.csv("BabayanEtAl_VirusData.csv",header=T)
fis<-read.csv("featureImportance_arthropodBorne.csv",header=T)

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)
gen.feats<-c(dinucs,cps,aa.codon.bias)
total.feats<-length(gen.feats)

f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector",gen.feats)]
f1$response<-factor(f1$Vector.borne)

# Remove viruses with unknown vector status
vecorphans<-subset(f1,f1$Vector.borne=="?")
vecorphans<-vecorphans[with(vecorphans,order(Viral.group,Genbank.accession,decreasing=T)),]
f_v<-subset(f1,f1$Vector.borne!="?")

# Choose which level of host group to analyze
f_v$response<-f_v$Vector.borne
f_v$response<-factor(f_v$response)

# Feature selection
nfeats<-25
totalfeats<-length(fis$vimean)
f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
gen.feats<-as.character(fis$X[f])
f_v2<-f_v[,c("Virus.name","response","Viral.group","Genbank.accession",gen.feats)] # simplify dataset to needed columns
ntax<-length(unique(f_v2$response))
bp<-c("X0","X1")
f_v2$Genbank.accession2<-do.call(rbind,strsplit(as.character(f_v2$Genbank.accession),"[.]"))[,1]
nfeatures<-length(gen.feats)

# Cleanup
rm(fis,f1,f_v)

# Train many models
set.seed(78910)
s<-.7 # proportion in the training set
nloops=600
accuracy.v<-c()
pc.accuracy<-matrix(nrow=nloops,ncol=2)
test.record<-matrix(nrow=80,ncol=nloops)
lr<-c()
md<-c()
sr<-c()
csr<-c()
nt<-c()
mr<-c()
vimps<-matrix(nrow=nfeatures,ncol=nloops)

for (i in 1:nloops){
  # 70:15:15 stratified split
    trains<-f_v2 %>% group_by(response) %>%
      filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(s*length(unique(Genbank.accession)))))
    testval<-subset(f_v2,!(f_v2$Genbank.accession %in% trains$Genbank.accession))
    vals<-testval %>% group_by(response) %>%
      filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
    tests<-subset(testval,!(testval$Genbank.accession %in% vals$Genbank.accession))
    ntest<-dim(tests)[1]

    trains<-droplevels(trains)
    tests<-droplevels(tests)
    vals<-droplevels(vals)
    test.record[,i]<-as.character(tests$Genbank.accession)

    # Build training, test, optimization and validation datasets
    set<-c("response",gen.feats)
    f1_train<-trains[,c(set)]

    testID<-tests$Genbank.accession
    f1_test<-tests[,c(set)]

    valID<-vals$Genbank.accession
    set<-c("response",gen.feats)
    f1_val<-vals[,c(set)]

    set<-c(gen.feats)
    f1_orphan<-vecorphans[,c(set)]


    # Convert to h2o data frames
    train<-as.h2o(f1_train)
    test<-as.h2o(f1_test)
    val<-as.h2o(f1_val)
    orp<-as.h2o(f1_orphan)

    # Clean up
    rm(f1_train,f1_test,f1_orphan,f1_val)

    # Identify the response column
    y <- "response"

    # Identify the predictor columns
    x <- setdiff(names(train), y)

    # Convert response to factor
    train[,y] <- as.factor(train[,y])
    test[,y] <- as.factor(test[,y])
    val[,y] <- as.factor(val[,y])

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

    gbm_grid4 <- h2o.grid("gbm", x = x, y = y,
    grid_id = "gbm_grid4",
    training_frame = train,
    validation_frame = val,
    seed = 1,
    hyper_params = gbm_params,
    search_criteria = search_criteria)

    gbm_gridperf <- h2o.getGrid(grid_id = "gbm_grid4",
    sort_by = "accuracy",
    decreasing = TRUE)

    # Grab the model_id for the top GBM model
    best_gbm_model_id <- gbm_gridperf@model_ids[[1]]
    best_gbm <- h2o.getModel(best_gbm_model_id)
    perf <- h2o.performance(best_gbm, test)

    # write model (if required)
    #h2o.saveModel(best_gbm,path =paste("h2o_ab_selGen_model",i))

    # record with best settings
    lr[i]<-as.numeric(gbm_gridperf@summary_table[1,2]) # learn_rate
    sr[i]<-as.numeric(gbm_gridperf@summary_table[1,6]) # sample_rate
    md[i]<-as.numeric(gbm_gridperf@summary_table[1,3]) # maxdepth
    csr[i]<-as.numeric(gbm_gridperf@summary_table[1,1]) # col_sample_rate
    mr[i]<-as.numeric(gbm_gridperf@summary_table[1,4]) # col_sample_rate
    nt[i]<-as.numeric(gbm_gridperf@summary_table[1,5]) # col_sample_rate

    # Print confusion matrix
    cm1<-h2o.confusionMatrix(perf)
    nclass<-length(unique(trains$response))
    cm2<-cm1[1:nclass,1:nclass]
    cm<-as.matrix(cm2)

    norm_cm<-cm/rowSums(cm)
    accuracy.v[i]=sum(diag(cm))/sum(cm)
    pc.accuracy[i,]<-t(diag(cm)/rowSums(cm))

    write.csv(norm_cm,file=paste("h2o_ab_selGen_CM",i,".csv"))

    # Retreive feature importance
    vi <- h2o.varimp(best_gbm)
    data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
    vimps[,i]<-data2[,4]

    # Orphan predictions
    orp.pred <- h2o.predict(best_gbm, orp)
    df<-orp.pred[,c(2:(ntax+1))]
    df2<-as.data.frame(df)
    row.names(df2)<-vecorphans$Virus.name

    write.csv(df2,file=paste("AB_Orphans",i,".csv",sep="_"))

    # Test set predictions
    test.pred<-h2o.predict(best_gbm,test[,2:length(names(test))]) # REMOVE host COLUMN
    df2<-as.data.frame(test.pred)
    row.names(df2)<-testID
    write.csv(df2,file=paste("ABTestPred",i,".csv",sep="_"))

    # Clean up
    h2o.rm("gbm_grid4") # remove any previous grid to prevent error when appending
    h2o.rm(best_gbm)
    rm(gbm_grid4,gbm_gridperf,cm1,cm2,df2)
}

accs<-data.frame(accuracy.v,pc.accuracy,lr,sr,md,csr,mr,nt)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(vimps,file="ArthropodBorne_SelGen25_FI.csv",row.names = T)
write.csv(accs,file="ArthropodBorne_SelGen25_out.csv",row.names=F)
write.csv(test.record,file="ArthropodBorne_SelGen25_TestSets.csv",row.names=F)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
res<-table(tests$response)
chanceAccurate<-round(sum(prob*res),digits=0)
tot<-sum(res)
nullAcc<-chanceAccurate/tot
print(nullAcc)
