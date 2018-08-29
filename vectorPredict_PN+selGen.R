# Babayan, Orton & Streicker 
# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes 
# Vector host prediction from selected genomic features and phylogenetic neighborhoods
# https://www.h2o.ai/products/h2o/ 

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o)
library(dplyr) 
library(reshape2)
library(ape)
library(seqinr)
library(matrixStats)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Start h2o JVM
localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# Read data from file
f1<-read.csv(file="BabayanEtAl_VirusData.csv",header=T)
allP<-read.fasta(file ="BabayanEtAl_sequences.fasta",  seqtype = "DNA", as.string = TRUE, seqonly = F, strip.desc = T)
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

# Organize vector database
mosq<-subset(f_v,f_v$Vector=="mosquito")
midg<-subset(f_v,f_v$Vector=="midge")
tick<-subset(f_v,f_v$Vector=="tick")
sand<-subset(f_v,f_v$Vector=="sandfly")

orpSeqs<-allP[c(which(names(allP) %in% vecorphans$Genbank.accession))] # pick sequences in the training set
write.fasta(orpSeqs, names(orpSeqs), file.out="vecorphanDB.fasta", open = "w", nbchar = 100, as.string = T)

# Feature selection
nfeats<-100
totalfeats<-length(fis$vimean)
f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
gen.feats<-as.character(fis[f,1])
bp<-as.character(sort(unique(f_v$Vector)))
total.feats<-length(gen.feats)+length(bp)

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

  trainSeqs<-allP[c(which(names(allP) %in% trains$Genbank.accession))] # pick sequences in the training set
  testSeqs<-allP[c(which(names(allP) %in% tests$Genbank.accession))] # pick sequences in the training set
  valSeqs<-allP[c(which(names(allP) %in% vals$Genbank.accession))] # pick sequences in the training set

  # Write sequence data to local directory
  write.fasta(testSeqs, names(testSeqs), file.out="testDB.fasta", open = "w", nbchar = 100, as.string = T)
  write.fasta(trainSeqs, names(trainSeqs), file.out="trainDB.fasta", open = "w", nbchar = 100, as.string = T)
  write.fasta(valSeqs, names(valSeqs), file.out="valDB.fasta", open = "w", nbchar = 100, as.string = T)

  # BLAST
  system("makeblastdb -in trainDB.fasta -dbtype nucl -parse_seqids -out allTrainingDB",intern=F)
  
  # Blast orphans against training
  system("blastn -db allTrainingDB -query vecorphanDB.fasta -out vecOrphanOut.out -num_threads 10 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F, wait=T)
  # Blast validation against training
  system("blastn -db allTrainingDB -query valDB.fasta -out valOut.out -num_threads 10 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F, wait=T)
  # Blast test against training
  system("blastn -db allTrainingDB -query testDB.fasta -out testOut.out -num_threads 10 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=T)
  # Blast training against the training set (take top 5 hits)
  system("blastn -db allTrainingDB -query trainDB.fasta -out trainOut.out -num_threads 10 -outfmt 10 -max_target_seqs=6 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=T)
  
  # Summarize blast hits from training set
  allBlast<-read.csv(file="trainOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
  nvir<-length(unique(allBlast$query.acc.))
  virnames<-unique(allBlast$query.acc.)
  ecutoff<-1E-3
  j=1
  d<-subset(allBlast,allBlast$query.acc.==virnames[j])
  d2<-subset(d,d$X..identity<100)
  d2<-subset(d2,d2$X.evalue<ecutoff)
  
  # Assign equal probability across all hosts if there is no good blast hit
  for (z in 1:1){
    if (nrow(d2)==0){
      blast.uc<-rep(1/ntax,ntax)
      blast.uc<-data.frame(t(blast.uc))
      colnames(blast.uc)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,blast.uc)}
    else {
      dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,hosts)}}
  
  for (j in 2:nvir){
    d<-subset(allBlast,allBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X..identity<100)
    d2<-subset(d2,d2$X.evalue<ecutoff)
    # Assign equal probability across all hosts if there is no good blast hit
    if (nrow(d2)==0){
      blast.uc.s<-rep(1/ntax,ntax)
      blast.uc.s<-data.frame(t(blast.uc.s))
      colnames(blast.uc.s)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc.s<-cbind(id,blast.uc.s)}
    else {
      dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(virnames[j])
      blast.uc.s<-cbind(id,hosts)}
    blast.uc<-rbind(blast.uc,blast.uc.s)}
  
  f1_train<-merge(trains,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
  set<-c("response",gen.feats,bp) 
  f1_train<-f1_train[,c(set)] 
  
  # Summarize blast hits from test set
  testBlast<-read.csv(file="testOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
  nvir<-length(unique(testBlast$query.acc.))
  virnames<-unique(testBlast$query.acc.)
  ecutoff<-1E-3
  j=1
  d<-subset(testBlast,testBlast$query.acc.==virnames[j])
  d2<-subset(d,d$X.evalue<ecutoff)
  
  # Assign equal probability across all hosts if there is no good blast hit
  for (z in 1:1){
    if (nrow(d2)==0){
      blast.uc<-rep(1/ntax,ntax)
      blast.uc<-data.frame(t(blast.uc))
      colnames(blast.uc)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,blast.uc)}
    else {
      dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,hosts)}}
  
  for (j in 2:nvir){
    d<-subset(testBlast,testBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)
    
    # Assign equal probability across all hosts if there is no good blast hit
    if (nrow(d2)==0){
      blast.uc.s<-rep(1/ntax,ntax)
      blast.uc.s<-data.frame(t(blast.uc.s))
      colnames(blast.uc.s)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc.s<-cbind(id,blast.uc.s) }
    else {
      dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(d$query.acc.[1])
      blast.uc.s<-cbind(id,hosts)}
    blast.uc<-rbind(blast.uc,blast.uc.s)}
  
  f1_test<-merge(tests,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
  testID<-f1_test$Virus.name
  set<-c("response",gen.feats,bp) 
  f1_test<-f1_test[,c(set)] 
  
  # Optimization set
  valBlast<-read.csv(file="valOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
  nvir<-length(unique(valBlast$query.acc.))
  virnames<-unique(valBlast$query.acc.)
  ecutoff<-1E-3
  j=1
  d<-subset(valBlast,valBlast$query.acc.==virnames[j])
  d2<-subset(d,d$X.evalue<ecutoff)
  
  # Assign equal probability across all hosts if there is no good blast hit
  for (z in 1:1){
    if (nrow(d2)==0){
      blast.uc<-rep(1/ntax,ntax)
      blast.uc<-data.frame(t(blast.uc))
      colnames(blast.uc)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,blast.uc)}
    else {
      dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,hosts)}}
  
  for (j in 2:nvir){
    d<-subset(valBlast,valBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)
    # Assign equal probability across all hosts if there is no good blast hit
    if (nrow(d2)==0){
      blast.uc.s<-rep(1/ntax,ntax)
      blast.uc.s<-data.frame(t(blast.uc.s))
      colnames(blast.uc.s)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc.s<-cbind(id,blast.uc.s) }
    else {
      dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(d$query.acc.[1])
      blast.uc.s<-cbind(id,hosts)}
    blast.uc<-rbind(blast.uc,blast.uc.s)}
  
  f1_val<-merge(vals,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
  valID<-f1_val$Virus.name
  set<-c("response",gen.feats,bp)
  f1_val<-f1_val[,c(set)] 
  
  # Vector orphans
  oBlast<-read.csv(file="vecOrphanOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
  nvir<-length(unique(oBlast$query.acc.))
  virnames<-unique(oBlast$query.acc.)
  ecutoff<-1E-3
  j=1
  d<-subset(oBlast,oBlast$query.acc.==virnames[j])
  d2<-subset(d,d$X.evalue<ecutoff)
  
  # Assign equal probability across all hosts if there is no good blast hit
  for (z in 1:1){
    if (nrow(d2)==0){
      blast.uc<-rep(1/ntax,ntax)
      blast.uc<-data.frame(t(blast.uc))
      colnames(blast.uc)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,blast.uc)}
    else {
      dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(virnames[j])
      blast.uc<-cbind(id,hosts)}}
  
  for (j in 2:nvir){
    d<-subset(oBlast,oBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)
    # Assign equal probability across all hosts if there is no good blast hit
    if (nrow(d2)==0){
      blast.uc.s<-rep(1/ntax,ntax)
      blast.uc.s<-data.frame(t(blast.uc.s))
      colnames(blast.uc.s)<-sort(unique(trains$response))
      id<-as.character(virnames[j])
      blast.uc.s<-cbind(id,blast.uc.s) }
    else {
      dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
      dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
      hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
      hosts[is.na(hosts)]<-0
      hosts<-t(data.frame(hosts))
      hosts<-data.frame(hosts)
      id<-as.character(d$query.acc.[1])
      blast.uc.s<-cbind(id,hosts)}
    blast.uc<-rbind(blast.uc,blast.uc.s)}
  
  f1_orphan<-merge(vecorphans,blast.uc,by.x="Genbank.accession",by.y="id",all.x=T,all.y=T,sort=F)
  set<-c(gen.feats,bp) 
  f1_orphan<-f1_orphan[,c(set)] 
  
  # Convert to h2o data frames    
  train<-as.h2o(f1_train)
  val<-as.h2o(f1_val)
  test<-as.h2o(f1_test)
  orpvec<-as.h2o(f1_orphan)

  # Clean up
  rm(f1_test,f1_train,f1_val,f1_orphan,allBlast,valBlast,testBlast,oBlast,trainSeqs,testSeqs,valSeqs)
  
  # Identity the response column
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
  
  # Write model (if required)
  #h2o.saveModel(best_gbm,path =paste("vectorPredict_PNSG_model",i))
  
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
  
  h2o.rm("gbm_grid4") 
  rm(gbm_grid4,gbm_gridperf,best_gbm,train,val,test,orpvec)
}

accs<-data.frame(accuracy.v,pc.accuracy,lr,sr,md,csr,mr,nt)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(accs,file="Vector_PNSelGen100_out.csv",row.names = F)
write.csv(test.record,file="Vector_PNSelGen100_TestSets.csv",row.names = F)
write.csv(vimps,file="Vector_PNSelGen100_FI.csv",row.names = T)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
vecs<-table(tests$response)
chanceAccurate<-round(sum(prob*vecs),digits=0)
tot<-sum(vecs)
nullAcc<-chanceAccurate/tot
print(nullAcc)