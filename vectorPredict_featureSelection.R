"""
Babayan, Orton & Streicker

Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

-- Genomic bias feature selection for vector taxon model
"""

rm(list=ls())

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o/
library(dplyr)
library(reshape2)
library(matrixStats)

setwd("") # Set local working directory where files are located

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
f1$response<-factor(f1$Vector)

# Select viral groups with vectors and remove features with little variation
f2<-subset(f1,f1$Viral.group=="Bunyavirus"|f1$Viral.group=="Flavivirus"|f1$Viral.group=="Rhabdovirus"|f1$Viral.group=="Togavirus")

f2_genomic<-f2[,7:4235]
f3<-f2_genomic[, sapply(f2_genomic, function(col) length(unique(col))) > 60]
ngenomic<-length(f3)
all.gen<-names(f3)
f4<-cbind(f2[1:6],f3)

# Select viruses with 4 vector classes
f5<-subset(f4,f4$Vector.borne==1)
f6<-subset(f5,f5$Vector=="mosquito"|f5$Vector=="midge"|f5$Vector=="sandfly"|f5$Vector=="tick")
f_v<-droplevels(f6)

# Clean up unneeded files
rm(f1,f2,f2_genomic,f3,f4,f5,f6)

# organize vector database
mosq<-subset(f_v,f_v$Vector=="mosquito")
midg<-subset(f_v,f_v$Vector=="midge")
tick<-subset(f_v,f_v$Vector=="tick")
sand<-subset(f_v,f_v$Vector=="sandfly")

s<-.7 # proportion in the training set

# Evaluate feature importance over 50 training sets (70% each with 5x cross-validation)
set.seed(78910)
nloops<-50
vimps<-matrix(nrow=ngenomic,ncol=nloops)

for (i in 1:nloops){
  mosq_sel<-mosq[sample(1:nrow(mosq),20),] # Downsample mosquito viruses
  f_all<-rbind(mosq_sel,midg,tick,sand)

  trains<-f_all %>% group_by(Vector) %>%
    filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(s*length(unique(Genbank.accession)))))

  set<-c("Vector",all.gen)
  f1_train<-trains[,c(set)]

  # Build h2o data frames
  train<-as.h2o(f1_train)

  # Identity the response column
  y <- "Vector"

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
                    keep_cross_validation_predictions=F)

  # Retreive feature importance
  vi <- h2o.varimp(model)
  data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
  vimps[,i]<-data2[,4] # percentage importance
  h2o.rm(model)
  rm(train,f1_train,vi,f_all,trains)
}

# Average feature importance across all training sets
row.names(vimps)<-data2$variable
vimean<-rowMeans2(vimps)
visd<-rowSds(vimps)
vistderr<-visd/sqrt(nloops)
vimps<-cbind(vimps,vimean,visd,vistderr)
vimps<- vimps[order(vimps[,nloops+1],decreasing=FALSE),]

# Write to file
write.csv(vimps,file="featureImportance_vector.csv",row.names = T)
