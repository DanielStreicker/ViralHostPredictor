# Babayan, Orton & Streicker 
# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes 
# Arthropod-bone transmission prediction from selected genomic features and phylogenetic neighborhoods
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
f_v2<-f_v[,c("Virus.name","response","Viral.group","Genbank.accession",gen.feats)] 
ntax<-length(unique(f_v2$response))
bp<-c("X0","X1")
f_v2$Genbank.accession2<-do.call(rbind,strsplit(as.character(f_v2$Genbank.accession),"[.]"))[,1]
nfeatures<-length(gen.feats)+ntax

# Read nucleotide sequences
allP<-read.fasta(file ="BabayanEtAl_sequences.fasta",  seqtype = "DNA", as.string = TRUE, seqonly = F, strip.desc = T)

# Write orphans
vecOrp<-allP[c(which(names(allP) %in% vecorphans$Genbank.accession))]
write.fasta(vecOrp,names(vecOrp),file.out="vecOrphanDB.fasta", open = "w", nbchar = 100, as.string = T)

# Clean up
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
    testval<-subset(f_v2,!(f_v2$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set
    vals<-testval %>% group_by(response) %>%
      filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
    tests<-subset(testval,!(testval$Genbank.accession %in% vals$Genbank.accession)) # ref numbers in testval set absent from test set    
    ntest<-dim(tests)[1]
  
    trains<-droplevels(trains)
    tests<-droplevels(tests)
    vals<-droplevels(vals)
    test.record[,i]<-as.character(tests$Genbank.accession)
  
    # Select sequences and write to file
    trainSeqs<-allP[c(which(names(allP) %in% trains$Genbank.accession))] # pick sequences in the training set
    testSeqs<-allP[c(which(names(allP) %in% tests$Genbank.accession))] # pick sequences in the test set
    valSeqs<-allP[c(which(names(allP) %in% vals$Genbank.accession))] # pick sequences in the optimization set
    
    write.fasta(testSeqs, names(testSeqs), file.out="testDB.fasta", open = "w", nbchar = 100, as.string = T)
    write.fasta(trainSeqs, names(trainSeqs), file.out="trainDB.fasta", open = "w", nbchar = 100, as.string = T)
    write.fasta(valSeqs, names(valSeqs), file.out="valDB.fasta", open = "w", nbchar = 100, as.string = T)
    
    # BLAST
    system("makeblastdb -in trainDB.fasta -dbtype nucl -parse_seqids -out allTrainingDB",intern=F)

    # Blast test against training
    system("blastn -db allTrainingDB -query testDB.fasta -out testOut.out -num_threads 8 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Blast validation against training
    system("blastn -db allTrainingDB -query valDB.fasta -out valOut.out -num_threads 8 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Blast orphan against training
    system("blastn -db allTrainingDB -query vecOrphanDB.fasta -out vecOrphanOut.out -num_threads 8 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)

    # Blast training against the training set (take top 5 hits)
    system("blastn -db allTrainingDB -query trainDB.fasta -out trainOut.out -num_threads 8 -outfmt 10 -max_target_seqs=6 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=TRUE)

    # Summarize blast hits from training set
    allBlast<-read.csv(file="trainOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(allBlast$query.acc.))
    virnames<-unique(allBlast$query.acc.)
    ecutoff<-1E-3
    j=1
    d<-subset(allBlast,allBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X..identity<100)
    d2<-subset(d2,d2$X.evalue<ecutoff)

    # Assign equal probability across all hosts if there is no acceptable blast hit
    for (z in 1:1){
        if (nrow(d2)==0){
            blast.uc<-rep(1/ntax,ntax)
            blast.uc<-data.frame(t(blast.uc))
            colnames(blast.uc)<-sort(unique(trains$response))
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,blast.uc)
            colnames(blast.uc)<-c("id","X0","X1")}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,hosts)
            colnames(blast.uc)<-c("id","X0","X1")}}

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
            blast.uc.s<-cbind(id,blast.uc.s)
            colnames(blast.uc.s)<-c("id","X0","X1")}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc.s<-cbind(id,hosts)
            colnames(blast.uc.s)<-c("id","X0","X1")}
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
            blast.uc<-cbind(id,blast.uc)
            colnames(blast.uc)<-c("id","X0","X1")}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,hosts)
            colnames(blast.uc)<-c("id","X0","X1")}}

    for (j in 2:nvir){
        d<-subset(testBlast,testBlast$query.acc.==virnames[j])
        d2<-subset(d,d$X.evalue<ecutoff)
        # Assign equal probability across all hosts if there is no good blast hit
        if (nrow(d2)==0){
            blast.uc.s<-rep(1/ntax,ntax)
            blast.uc.s<-data.frame(t(blast.uc.s))
            colnames(blast.uc.s)<-sort(unique(trains$response))
            id<-as.character(virnames[j])
            blast.uc.s<-cbind(id,blast.uc.s)
            colnames(blast.uc.s)<-c("id","X0","X1") }
        else {
            dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(d$query.acc.[1])
            blast.uc.s<-cbind(id,hosts)
            colnames(blast.uc.s)<-c("id","X0","X1")}
        blast.uc<-rbind(blast.uc,blast.uc.s)}

    f1_test<-merge(tests,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
    testID<-f1_test$Genbank.accession
    set<-c("response",gen.feats,bp) 
    f1_test<-f1_test[,c(set)] 

    # Summarize blast hits from optimization set
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
        blast.uc<-cbind(id,blast.uc)
        colnames(blast.uc)<-c("id","X0","X1")}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,hosts)
        colnames(blast.uc)<-c("id","X0","X1")}}
    
    for (j in 2:nvir){
      d<-subset(valBlast,valBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      # Assign equal probability across all hosts if there is no good blast hit
      if (nrow(d2)==0){
        blast.uc.s<-rep(1/ntax,ntax)
        blast.uc.s<-data.frame(t(blast.uc.s))
        colnames(blast.uc.s)<-sort(unique(trains$response))
        id<-as.character(virnames[j])
        blast.uc.s<-cbind(id,blast.uc.s)
        colnames(blast.uc.s)<-c("id","X0","X1") }
      else {
        dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(d$query.acc.[1])
        blast.uc.s<-cbind(id,hosts)
        colnames(blast.uc.s)<-c("id","X0","X1")}
      blast.uc<-rbind(blast.uc,blast.uc.s)}
    
    f1_val<-merge(vals,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
    valID<-f1_val$Genbank.accession
    set<-c("response",gen.feats,bp)
    f1_val<-f1_val[,c(set)] 
    
    # Summarize blast hits from orphan set
    oBlast<-read.csv(file="vecOrphanOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(oBlast$query.acc.))
    virnames<-unique(oBlast$query.acc.)
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
        blast.uc<-cbind(id,blast.uc)
        colnames(blast.uc)<-c("id","X0","X1")}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=F)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,hosts)
        colnames(blast.uc)<-c("id","X0","X1")}}

    for (j in 2:nvir){
      d<-subset(oBlast,oBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      # Assign equal probability across all hosts if there is no good blast hit
      if (nrow(d2)==0){
        blast.uc.s<-rep(1/ntax,ntax)
        blast.uc.s<-data.frame(t(blast.uc.s))
        colnames(blast.uc.s)<-sort(unique(trains$response))
        id<-as.character(virnames[j])
        blast.uc.s<-cbind(id,blast.uc.s)
        colnames(blast.uc.s)<-c("id","X0","X1") }
      else {
        dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$response,sum,na.rm=T)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(d$query.acc.[1])
        blast.uc.s<-cbind(id,hosts)
        colnames(blast.uc.s)<-c("id","X0","X1")}
      blast.uc<-rbind(blast.uc,blast.uc.s)}

    f1_orphan<-merge(vecorphans,blast.uc,by.x="Genbank.accession",by.y="id",all.x=T,all.y=T,sort=F)
    names(f1_orphan)[c(4237,4238)]<-c("X0","X1")
    set<-c(gen.feats,"X0","X1")
    f1_orphan<-f1_orphan[,c(set)] 

    # Convert to h2o data frames    
    train<-as.h2o(f1_train)
    test<-as.h2o(f1_test)
    val<-as.h2o(f1_val)
    orp<-as.h2o(f1_orphan)
    
    # Clean up
    rm(f1_orphan,f1_val,f1_test,f1_train,trainSeqs,testSeqs,valSeqs)
    
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
    
    # Write model (if required)
    #h2o.saveModel(best_gbm,path =paste("h2o_ab_selGen.blastn_model",i))
    
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
    
    write.csv(norm_cm,file=paste("h2o_ab_selGen.blastn_CM",i,".csv"))
    
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
    
    # Test set 
    test.pred<-h2o.predict(best_gbm,test[,2:length(names(test))]) # REMOVE host COLUMN
    df2<-as.data.frame(test.pred)
    row.names(df2)<-testID
    write.csv(df2,file=paste("ABTestPred",i,".csv",sep="_"))  
    
    h2o.rm(c("gbm_grid4",best_gbm))
    rm(vals,trainSeqs,testSeqs,valSeqs,oBlast,valBlast,testBlast,allBlast,f1_train,f1_test,f1_val,f1_orphan)
}

accs<-data.frame(accuracy.v,pc.accuracy,lr,sr,md,csr,mr,nt)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(vimps,file="ArthropodBorne_PNSelGen25_FI.csv",row.names = T)
write.csv(accs,file="ArthropodBorne_PNSelGen25_out.csv",row.names=F)
write.csv(test.record,file="ArthropodBorne_PNSelGen25_TestSets.csv",row.names=F)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
res<-table(tests$response)
chanceAccurate<-round(sum(prob*res),digits=0)
tot<-sum(res)
nullAcc<-chanceAccurate/tot
print(nullAcc)