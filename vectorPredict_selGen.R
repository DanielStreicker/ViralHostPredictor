"""
Babayan, Orton & Streicker

Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

-- Vector host prediction from selected genomic features
"""

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o/
library(dplyr)
library(reshape2)
library(matrixStats)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Start h2o JVM
localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# Read data from file
f1<-read.csv(file="BabayanEtAl_VirusData.csv",header=T)
fis<-read.csv(file="featureImportance_vector.csv",header=T)

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)

# Remove orphans
f2<-subset(f1,f1$Vector.borne==1)
f3<-subset(f2,f2$Vector=="mosquito"|f2$Vector=="midge"|f2$Vector=="sandfly"|f2$Vector=="tick")
f_v<-droplevels(f3)
f_v$response<-f_v$Vector
vecorphans<-subset(f1,f1$Genbank.accession=="KF186496.1"|f1$Genbank.accession=="NC_008718.1"|f1$Genbank.accession=="NC_005039.1"|f1$Genbank.accession=="NC_026624.1"|f1$Genbank.accession=="NC_022755.1"|f1$Genbank.accession=="JX297815.1"|f1$Genbank.accession=="NC_025341.1"|f1$Genbank.accession=="NC_007020.1"|f1$Genbank.accession=="NC_025391.1"|f1$Genbank.accession=="NC_026623.1"|f1$Genbank.accession=="NC_002526.1"|f1$Genbank.accession=="NC_009026.2"|f1$Genbank.accession=="NC_025396.1"|f1$Genbank.accession=="NC_025358.1"|f1$Genbank.accession=="NC_025253.1"|f1$Genbank.accession=="NC_023812.1"|f1$Genbank.accession=="NC_015375.1"|f1$Genbank.accession=="NC_007020.1"|f1$Genbank.accession=="KM408491.1")
vecorphans<-vecorphans[with(vecorphans,order(Viral.group,Genbank.accession,decreasing=T)),]

# organize vector database
mosq<-subset(f_v,f_v$Vector=="mosquito")
midg<-subset(f_v,f_v$Vector=="midge")
tick<-subset(f_v,f_v$Vector=="tick")
sand<-subset(f_v,f_v$Vector=="sandfly")

# feature selection
nfeats<-100
totalfeats<-length(fis$vimean)
f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
gen.feats<-as.character(fis[f,1])
bp<-as.character(sort(unique(f_v$Vector)))
total.feats<-length(gen.feats)

# Clean up
rm(f1,f2,f3)

# Train many models
set.seed(78910)
s<-.6 # Proportion in the training set
nloops<-250
lr<-c()
md<-c()
sr<-c()
csr<-c()
nt<-c()
mr<-c()
accuracy.v<-c()
pc.accuracy<-matrix(nrow=nloops,ncol=4)
test.record<-matrix(nrow=46,ncol=nloops)

vimps<-matrix(nrow=total.feats,ncol=nloops)
ntax=length(bp)

for (i in 1:nloops){
  mosq_sel<-mosq[sample(1:nrow(mosq),20),]   # downsample mosquito viruses
  f_all<-rbind(mosq_sel,midg,tick,sand)

  trains<-f_all %>% group_by(response) %>%
    filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(s*length(unique(Genbank.accession)))))
  testval<-subset(f_v,!(f_v$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set

  vals<-testval %>% group_by(response) %>%
    filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.3*length(unique(Genbank.accession)))))
  tests<-subset(testval,!(testval$Genbank.accession %in% vals$Genbank.accession)) # ref numbers in testval set absent from test set

  trains<-droplevels(trains)
  tests<-droplevels(tests)
  vals<-droplevels(vals)

  ntest<-dim(tests)[1]
  test.record[,i]<-as.character(tests$Genbank.accession)

  # Training set
  set<-c("response",gen.feats)
  f1_train<-trains[,c(set)]

  # Test set
  testID<-tests$Virus.name
  f1_test<-tests[,c(set)]

  # Optimization set
  valID<-vals$Virus.name
  f1_val<-vals[,c(set)]

  # Vector orphans
  set<-c(gen.feats)
  f1_orphan<-vecorphans[,c(set)]

  # Convert to h2o data frames
  train<-as.h2o(f1_train)
  val<-as.h2o(f1_val)
  test<-as.h2o(f1_test)
  orpvec<-as.h2o(f1_orphan)

  rm(f1_test,f1_train,f1_val,f1_orphan)

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
                        max_models = 500)

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

  # Record best settings
  lr[i]<-as.numeric(gbm_gridperf@summary_table[1,2]) # learn_rate
  sr[i]<-as.numeric(gbm_gridperf@summary_table[1,6]) # sample_rate
  md[i]<-as.numeric(gbm_gridperf@summary_table[1,3]) # maxdepth
  csr[i]<-as.numeric(gbm_gridperf@summary_table[1,1]) # col_sample_rate
  mr[i]<-as.numeric(gbm_gridperf@summary_table[1,4]) # min rows
  nt[i]<-as.numeric(gbm_gridperf@summary_table[1,5]) # n trees

  # Print confusion matrix
  cm1<-h2o.confusionMatrix(perf)
  nclass<-length(unique(trains$response))
  cm2<-cm1[1:nclass,1:nclass]
  cm<-as.matrix(cm2)

  norm_cm<-cm/rowSums(cm)
  accuracy.v[i]=sum(diag(cm))/sum(cm)
  pc.accuracy[i,]<-t(diag(cm)/rowSums(cm))

  rownames(norm_cm)<-c("Midge","Mosquito","Sandfly","Tick")
  colnames(norm_cm)<-c("Midge","Mosquito","Sandfly","Tick")
  write.csv(norm_cm,file=paste("h2o_Vector_dinuc.blastn_CM_",i,".csv"))

  # Retreive feature importance
  vi <- h2o.varimp(best_gbm)
  data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
  vimps[,i]<-data2[,4]

  # Orphan predictions
  orp.pred <- h2o.predict(best_gbm, orpvec)
  df<-orp.pred[,c(2:5)]
  df2<-as.data.frame(df)
  row.names(df2)<-vecorphans$Virus.name
  colnames(df2)<-c("Midge","Mosquito","Sandfly","Tick")

  write.csv(df2,file=paste("VectorOrphans",i,".csv",sep="_"))

  # Test set predictions
  test.pred<-h2o.predict(best_gbm,test[,2:length(names(test))]) # REMOVE host COLUMN
  df2<-as.data.frame(test.pred)
  row.names(df2)<-testID
  write.csv(df2,file=paste("VectorTestPred",i,".csv",sep="_"))

  # Clean up
  h2o.rm("gbm_grid4")
  rm(gbm_grid4,gbm_gridperf,best_gbm,train,val,test,orpvec,test.pred,orp.pred,vi)
}

accs<-data.frame(accuracy.v,pc.accuracy,lr,sr,md,csr,mr,nt)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(accs,file="Vector_SelGen100_out.csv",row.names = F)
write.csv(test.record,file="Vector_SelGen100_TestSets.csv",row.names = F)
write.csv(vimps,file="Vector_SelGen100_FI.csv",row.names = T)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
vecs<-table(tests$response)
chanceAccurate<-round(sum(prob*vecs),digits=0)
tot<-sum(vecs)
nullAcc<-chanceAccurate/tot
print(nullAcc)