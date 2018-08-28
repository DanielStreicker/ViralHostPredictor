# Babayan, Orton & Streicker 
# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes 
# Genomic bias feature selection for arthropod-borne transmission model
# https://www.h2o.ai/products/h2o/ 

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o)
library(dplyr)
library(reshape2)
library(matrixStats)

localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# Data
f1<-read.csv("BabayanEtAl_VirusData.csv")

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)
gen.feats<-c(dinucs,cps,aa.codon.bias)
total.feats<-length(gen.feats)

f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector",gen.feats)] 
f1$response<-factor(f1$Vector.borne)

# Remove viruses with unknown arthropod-borne status
f2<-subset(f1,f1$Vector.borne!="?")
f_v<-droplevels(f2)

# Clean up unneeded files
rm(f1,f2)

# Evaluate patterns over many training sets
set.seed(78910)
s<-.7 # Proportion in the training set
nloops<-50
vimps<-matrix(nrow=total.feats-2,ncol=nloops)

for (i in 1:nloops){
  trains<-f_v %>% group_by(Vector.borne) %>%
    filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(s*length(unique(Genbank.accession)))))

  set<-c("Vector.borne",gen.feats)
  f1_train<-trains[,c(set)]

  # Build h2o data frames
  train<-as.h2o(f1_train)

  # Identity the response column
  y <- "Vector.borne"

  # Identify the predictor columns
  x <- setdiff(names(train), y)

  # Convert response to factor
  train[,y] <- as.factor(train[,y])
  
  # GBM with 5x cross validation of training set, test set is not used
  model <- h2o.gbm(x = x,
                    y = y,
                    training_frame = train,
                    ntrees = 150,
                    learn_rate = .1,
                    sample_rate = 1,
                    max_depth = 10,
                    col_sample_rate_per_tree = 1,
                    seed = 123,
                    nfolds = 5,
                    keep_cross_validation_predictions=T)
  
  vi <- h2o.varimp(model)
  data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
  vimps[,i]<-data2[,4] # "percentage" importance
  h2o.rm(model)
  rm(train,f1_train,vi)
}

row.names(vimps)<-data2$variable

# Average feature importance across all training sets
vimean<-rowMeans2(vimps)
visd<-rowSds(vimps)
vistderr<-visd/sqrt(nloops)
vimps<-cbind(vimps,vimean,visd,vistderr)
vimps<- vimps[order(vimps[,nloops+1],decreasing=FALSE),]

# Write to file 
write.csv(vimps,file="featureImportance_arthropodBorne.csv",row.names = T)