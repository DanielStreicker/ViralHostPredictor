# This script analyzes the output of h2o.gbm test set predictions across many models

rm(list=ls())

library(matrixStats)
library(ggplot2)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# choose file set to use
#ttype<-"vector"
#ttype<-"reservoir"
ttype<-"vector.borne"

# choose working directory
ntest=80
wd1<-"" ## working directory name where results are stored
ptype<-as.character("VBTestPred")
setwd(wd1)

myfiles<-c()
nmodels<-600 # the number of files in the working directory to analyze
for(i in (1:nmodels)){
  myfiles[i]<-paste(ptype,i,".csv",sep="_")
}
testPredFiles = do.call(rbind, lapply(myfiles, function(x) read.csv(x, stringsAsFactors = FALSE)))

## assign overall model accuracy to each prediction
accs<-read.csv(file="",header=T) # validation set accuracies from nmodels splits of the data
acc.vector<-rep(accs$accuracy.v[1],ntest) # the accuracy of each model
for (i in 2:nmodels){
  av2<-rep(accs$accuracy.v[i],ntest)  
  acc.vector<-c(acc.vector,av2)
}
testPredFiles$acc<-acc.vector

allP.vb.psg<-testPredFiles
min(table(allP.vb.psg$X)) ## lowest number of observations of any virus in test set
range(table(allP.vb.psg$X)) ## range of number of observations of any virus in test set
median(table(allP.vb.psg$X)) ## median number of observations of any virus in test set

# generate model averaged predictions for class for each virus
t.model<-.75 ## threshold on predictions to use 1-t.model = the proportion of predictions used

## PREDICTIONS USING COMBINED PN AND GENOMIC MODEL
if (ttype=="vector.borne"){ allP_bsg<-allP.vb.psg}
nclass<-ncol(allP_bsg)-3
nvirus<-length(unique(allP_bsg$X))
vnames<-unique(allP_bsg$X)
vir.m<-matrix(nrow=nvirus,ncol=nclass)
vir.sd<-matrix(nrow=nvirus,ncol=nclass)
nc2<-nclass+2

for(i in 1:nvirus){
  vir<-subset(allP_bsg,allP_bsg$X==vnames[i])
  ## pick top 25% of models
  vir<-subset(vir,vir$acc>=quantile(vir$acc,t.model))
  vir.m[i,1:nclass]<-colMeans(vir[3:nc2])
  vir.sd[i,1:nclass]<-colSds(as.matrix(vir[3:nc2]))
}
row.names(vir.m)<-vnames
colnames(vir.m)<-names(vir[3:nc2])
row.names(vir.sd)<-vnames
colnames(vir.sd)<-names(vir[3:nc2])
vir.m.pnsg<-vir.m
vir.sd.pnsg<-vir.sd

# find the top prediction for each virus as the class with the highest support
vir.df<-data.frame(vir.m)
vir.df$pred<-colnames(vir.df)[max.col(vir.df,ties.method="first")]
vir.df$pp<-rowMaxs(as.matrix(vir.df[,c(1:nclass)])) # largest probability
vir.df$GenbankID<-rownames(vir.df)

if(ttype=="vector.borne"){vir.df.bsg.vb<-vir.df}

# merge predicted host classes with true host data
# load true data
f1<-read.csv(file="BabayanEtAl_VirusData.csv",header=T)

## feature definition
f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector")]

if (ttype=="vector.borne"){
f1$response<-factor(f1$Vector.borne)}

# evaluate model averaged accuracy 
## vector borne
mer.bsg.vb<-merge(vir.df.bsg.vb,f1,by.x="GenbankID",by.y="Genbank.accession",all.x = T,all.y=F)
mer.bsg.vb$pred2<-ifelse(mer.bsg.vb$pred=="p0",0,1)
mer.bsg.vb$score<-ifelse(mer.bsg.vb$pred2==mer.bsg.vb$response,1,0)
table(mer.bsg.vb$score)[2]/sum(table(mer.bsg.vb$score))

# ## WRITE OUT RESULTS FOR VECTOR BORNE
mer<-mer.bsg.vb
tab.sg.vb<-mer[,c(6,2:3,8,9,12,13)]
colnames(tab.sg.vb)<-c("Virus Name","No vector","Vector-borne","Viral group","True status","Predicted status","Score")
write.table(tab.sg.vb,file="VectorBorne_PNSG_StudyWideTestSet.csv",sep=",",row.names = F)

## per class accuracy
# arthropod borne
pcacc<-table(mer.bsg.vb$response,mer.bsg.vb$score)
perc<-pcacc[,2]/(pcacc[,1]+pcacc[,2])
perc<-na.exclude(perc)
nvirs<-table(mer.bsg$response)