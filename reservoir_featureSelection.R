"""
Babayan, Orton & Streicker

Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

-- Genomic bias feature selection for reservoir host model
"""

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o/
library(dplyr)
library(reshape2)
library(matrixStats)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# Data
f1<-read.csv(file="BabayanEtAl_VirusData.csv",header=T)

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)

f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector",dinucs,cps,aa.codon.bias)]
f1$response<-factor(f1$Reservoir)

# Remove orphans
f2<-subset(f1,f1$response!="Orphan")
f<-droplevels(f2)

# Group selection based on sample size thresholds
t<-15 # Minimum sample size per group
s<-.7 # Proportion in the training set
host.counts<-table(f$response)
min.t<-host.counts[host.counts>=t] # minimum number of viruses per host group
f_st3<-f[f$response %in% c(names(min.t)),]
f_st3<-droplevels(f_st3)

rm(f,f2,f1)

# Evaluate patterns over many training sets
set.seed(78910)
nloops<-50

vimps<-matrix(nrow=length(c(cps,dinucs,aa.codon.bias))-2,ncol=nloops) # two are removed because ATG.Bias and TGG.Bias are constant

for (i in 1:nloops){
    trains<-f_st3 %>% group_by(response) %>%
    filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
    trains<-droplevels(trains)

    set<-c("response",dinucs,cps,aa.codon.bias)
    f1_train<-trains[,c(set)]
    train<-as.h2o(f1_train)

    # Identity the response column
    y <- "response"

    # Identify the predictor columns
    x <- setdiff(names(train), y)

    # Convert response to factor
    train[,y] <- as.factor(train[,y])

    # GBM with 5x cross validation of training set, test set is not used
    model1 <- h2o.gbm(x = x,
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

    # Retreive feature importance
    vi <- h2o.varimp(model1)
    data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
    vimps[,i]<-data2[,4] # "percentage" importance
    h2o.rm(model1)
    rm(trains,train,f1_train,vi)
}

# Average feature importance across all training sets
row.names(vimps)<-data2$variable
vimean<-rowMeans2(vimps)
visd<-rowSds(vimps)
vimps<-cbind(vimps,vimean,visd)
vistderr<-visd/sqrt(nloops)
vimps<-cbind(vimps,vimean,visd,vistderr)
vimps<- vimps[order(vimps[,nloops+1],decreasing=FALSE),] # sort by mean feature importance

# Write to file
write.csv(vimps,file="featureImportance_reservoir.csv",row.names = T)
