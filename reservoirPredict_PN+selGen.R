# Babayan, Orton & Streicker 
# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes 
# Reservoir host prediction from selected genomic features and phylogenetic neighborhoods
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

# sample split of training/test to get counts in each 
trains<-f_st3 %>% group_by(Reservoir) %>%
  filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
testval<-subset(f_st3,!(f_st3$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set

optims<-testval %>% group_by(Reservoir) %>%
  filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
tests<-subset(testval,!(testval$Genbank.accession %in% optims$Genbank.accession)) # ref numbers in testval set absent from test set    
ntest<-dim(tests)[1]

# write orphan sequences
orp<-allP[c(which(names(allP) %in% orphans$Genbank.accession))]
write.fasta(orp,names(orp),file.out="orphanDB.fasta", open = "w", nbchar = 100, as.string = T)

# write rare host sequences
rar<-allP[c(which(names(allP) %in% rare$Genbank.accession))]
write.fasta(rar,names(rar),file.out="rareDB.fasta", open = "w", nbchar = 100, as.string = T)

# Remove unneeded files
rm(f,f1,f2,fis,rar,orp)


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
nfeatures<-length(gen.feats)+ntax
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
    
    # Select and write sequences to local directory
    trainSeqs<-allP[c(which(names(allP) %in% trains$Genbank.accession))] # pick sequences in the training set
    testSeqs<-allP[c(which(names(allP) %in% tests$Genbank.accession))] # pick sequences in the validation set
    optSeqs<-allP[c(which(names(allP) %in% optims$Genbank.accession))] # pick sequences in the optimization set
    write.fasta(testSeqs, names(testSeqs), file.out="testDB.fasta", open = "w", nbchar = 100, as.string = T)
    write.fasta(trainSeqs, names(trainSeqs), file.out="trainDB.fasta", open = "w", nbchar = 100, as.string = T)
    write.fasta(optSeqs, names(optSeqs), file.out="optDB.fasta", open = "w", nbchar = 100, as.string = T)
    
    # BLAST
    system("makeblastdb -in trainDB.fasta -dbtype nucl -parse_seqids -out allTrainingDB",intern=F)

    # Blast test against training
    system("blastn -db allTrainingDB -query testDB.fasta -out testOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Blast validation against training
    system("blastn -db allTrainingDB -query optDB.fasta -out optOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Orphan blast (take top 5 hits)
    system("blastn -db allTrainingDB -query orphanDB.fasta -out orphanOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Rare blast (take top 5 hits)
    system("blastn -db allTrainingDB -query rareDB.fasta -out rareOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    
    # Blast training against the training set (take top 5 hits)
    system("blastn -db allTrainingDB -query trainDB.fasta -out trainOut.out -num_threads 4 -outfmt 10 -max_target_seqs=6 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=TRUE)
    
    # Read in Blast output for training set
    allBlast<-read.csv(file="trainOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)

    # Summarize blast hits
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
            colnames(blast.uc)<-sort(unique(trains$Reservoir))
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,blast.uc)}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,hosts)}}

    for (j in 2:nvir){
        d<-subset(allBlast,allBlast$query.acc.==virnames[j])
        d2<-subset(d,d$X..identity<100)
        d2<-subset(d2,d2$X.evalue<ecutoff)
        if (nrow(d2)==0){
            blast.uc.s<-rep(1/ntax,ntax)
            blast.uc.s<-data.frame(t(blast.uc.s))
            colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
            id<-as.character(virnames[j])
            blast.uc.s<-cbind(id,blast.uc.s)}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc.s<-cbind(id,hosts)}
        blast.uc<-rbind(blast.uc,blast.uc.s)}

    f1_train<-merge(trains,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
    set<-c("Reservoir",gen.feats,bp) 
    f1_train<-f1_train[,c(set)] # this is the full training dataset with genomic features and blast probabilities

    # Summarize blast hits from test set
    testBlast<-read.csv(file="testOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(testBlast$query.acc.))
    virnames<-unique(testBlast$query.acc.)
    ecutoff<-1E-3
    j=1
    d<-subset(testBlast,testBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)

    for (z in 1:1){
        if (nrow(d2)==0){
            blast.uc<-rep(1/ntax,ntax)
            blast.uc<-data.frame(t(blast.uc))
            colnames(blast.uc)<-sort(unique(trains$Reservoir))
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,blast.uc)}
        else {
            dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(virnames[j])
            blast.uc<-cbind(id,hosts)}}

    for (j in 2:nvir){
        d<-subset(testBlast,testBlast$query.acc.==virnames[j])
        d2<-subset(d,d$X.evalue<ecutoff)
        if (nrow(d2)==0){
            blast.uc.s<-rep(1/ntax,ntax)
            blast.uc.s<-data.frame(t(blast.uc.s))
            colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
            id<-as.character(virnames[j])
            blast.uc.s<-cbind(id,blast.uc.s) }
        else {
            dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
            dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
            hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
            hosts[is.na(hosts)]<-0
            hosts<-t(data.frame(hosts))
            hosts<-data.frame(hosts)
            id<-as.character(d$query.acc.[1])
            blast.uc.s<-cbind(id,hosts)}
        blast.uc<-rbind(blast.uc,blast.uc.s)}

    f1_test<-merge(tests,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
    testID<-f1_test$Virus.name
    set<-c("Reservoir",gen.feats,bp)
    f1_test<-f1_test[,c(set)] 

    # Summarize blast hits from optimization set
    optBlast<-read.csv(file="optOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(optBlast$query.acc.))
    virnames<-unique(optBlast$query.acc.)
    ecutoff<-1E-3
    j=1
    d<-subset(optBlast,optBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)
    
    for (z in 1:1){
      if (nrow(d2)==0){
        blast.uc<-rep(1/ntax,ntax)
        blast.uc<-data.frame(t(blast.uc))
        colnames(blast.uc)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,blast.uc)}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,hosts)}}
    
    for (j in 2:nvir){
      d<-subset(optBlast,optBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      if (nrow(d2)==0){
        blast.uc.s<-rep(1/ntax,ntax)
        blast.uc.s<-data.frame(t(blast.uc.s))
        colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc.s<-cbind(id,blast.uc.s) }
      else {
        dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(d$query.acc.[1])
        blast.uc.s<-cbind(id,hosts)}
      blast.uc<-rbind(blast.uc,blast.uc.s)}
    
    f1_opt<-merge(optims,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
    optID<-f1_opt$Virus.name
    set<-c("Reservoir",gen.feats,bp)
    f1_opt<-f1_opt[,c(set)] 
    
    # Summarize blast hits from orphan set
    oBlast<-read.csv(file="orphanOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(oBlast$query.acc.))
    virnames<-unique(oBlast$query.acc.)
    j=1
    d<-subset(oBlast,oBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)

    for (z in 1:1){
      if (nrow(d2)==0){
        blast.uc<-rep(1/ntax,ntax)
        blast.uc<-data.frame(t(blast.uc))
        colnames(blast.uc)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,blast.uc)}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,hosts)}}

    for (j in 2:nvir){
      d<-subset(oBlast,oBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      if (nrow(d2)==0){
        blast.uc.s<-rep(1/ntax,ntax)
        blast.uc.s<-data.frame(t(blast.uc.s))
        colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc.s<-cbind(id,blast.uc.s) }
      else {
        dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(d$query.acc.[1])
        blast.uc.s<-cbind(id,hosts)}
      blast.uc<-rbind(blast.uc,blast.uc.s)}

    f1_orphan<-merge(orphans,blast.uc,by.x="Genbank.accession",by.y="id",all.x=T,all.y=T,sort=F)
    set<-c(gen.feats,bp) 
    f1_orphan<-f1_orphan[,c(set)] 

    # Summarize blast hits from rare virus set
    rBlast<-read.csv(file="rareOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(rBlast$query.acc.))
    virnames<-unique(rBlast$query.acc.)
    j=1
    d<-subset(rBlast,rBlast$query.acc.==virnames[j])
    d2<-subset(d,d$X.evalue<ecutoff)
    
    for (z in 1:1){
      if (nrow(d2)==0){
        blast.uc<-rep(1/ntax,ntax)
        blast.uc<-data.frame(t(blast.uc))
        colnames(blast.uc)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,blast.uc)}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j])
        blast.uc<-cbind(id,hosts)}}
    
    for (j in 2:nvir){
      d<-subset(rBlast,rBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      if (nrow(d2)==0){
        blast.uc.s<-rep(1/ntax,ntax)
        blast.uc.s<-data.frame(t(blast.uc.s))
        colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
        id<-as.character(virnames[j])
        blast.uc.s<-cbind(id,blast.uc.s) }
      else {
        dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
        dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(d$query.acc.[1])
        blast.uc.s<-cbind(id,hosts)}
      blast.uc<-rbind(blast.uc,blast.uc.s)}
    
    f1_rare<-merge(rare,blast.uc,by.x="Genbank.accession",by.y="id",all.x=T,all.y=T,sort=F)
    set<-c(gen.feats,bp) 
    f1_rare<-f1_rare[,c(set)] 
    
    
    # Convert to h2o data frames    
    train<-as.h2o(f1_train)
    test<-as.h2o(f1_test)
    opt<-as.h2o(f1_opt)
    orp<-as.h2o(f1_orphan)
    rar<-as.h2o(f1_rare)
    
    # Clean up
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
    
    # Print confusion matrix
    cm1<-h2o.confusionMatrix(perf)
    nclass<-length(unique(trains$Reservoir))
    cm2<-cm1[1:nclass,1:nclass]
    cm<-as.matrix(cm2)

    norm_cm<-cm/rowSums(cm)
    accuracy.st3[i]=sum(diag(cm))/sum(cm)
    pc.accuracy[i,]<-t(diag(cm)/rowSums(cm))
    
    write.csv(norm_cm,file=paste("Reservoir_selGen+PN_CM",i,".csv"))
    
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
    rm(oBlast,testBlast,allBlast,rBlast,trainSeqs,testSeqs,optSeqs,gbm_grid,best_gbm,train,test,opt,df2,optims)
}

accs<-data.frame(accuracy.st3,pc.accuracy,lr,sr,md,csr)
colnames(accs)[2:(ntax+1)]<-row.names(cm)
row.names(vimps)<-data2$variable

# Write results summaries
write.csv(vimps,file="Reservoir_PN+SelGen50_FI.csv",row.names = T)
write.csv(accs,file="Reservoir_PN+SelGen50_out.csv",row.names=F)

# Null model accuracy
prob<-table(trains$response)/sum(table(trains$response))
vecs<-table(tests$response)
chanceAccurate<-round(sum(prob*vecs),digits=0)
tot<-sum(vecs)
nullAcc<-chanceAccurate/tot
print(nullAcc)