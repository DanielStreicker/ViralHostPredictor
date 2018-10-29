# Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

_Simon A. Babayan, Richard J. Orton and Daniel G. Streicker_

## Background
A series of scripts and datasets described in Babayan _et al._ (2018) _Science_ which predict the reservoir hosts, existence of arthropod vectors and identity of arthropod vectors using gradient boosting machines implemented through h2o using R.

## File descriptions

### Datasets:

_BabayanEtAl_sequences.fasta_ contains coding sequences for all viruses used in the analyses

_EbolaTimeSeriesData.csv_ contains epidemiological data and genomic features for Zaire ebolaviruses sampled during the 2014-2016 West African outbreak

_BabayanEtAl_VirusData.csv_ contains reservoir host, arthropod-borne transmission status and vector taxa for all ssRNA viruses analyzed and features extracted from the genome of each virus

### R scripts:

_arthropodBorne_featureSelection.R_ Uses gradient boosting machines in h2o to estimate average feature importances for predicting arthropod-borne transmission across different training sets 

_arthropodBorne_PN+selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the arthropod-borne transmission status of each virus using phylogenetic neighborhoods and genomic features selected by _arthropodBorne_featureSelection.R_

_arthropodBorne_PN.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the arthropod-borne transmission status of each virus using phylogenetic neighborhoods

_arthropodBorne_selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the arthropod-borne transmission status of each virus using genomic features selected by _arthropodBorne_featureSelection.R_

_reservoir_featureSelection.R_ Uses gradient boosting machines in h2o to estimate average feature importances for predicting reservoir hosts across different training sets 

_reservoirPredict_PN+selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the reservoir host of each virus using phylogenetic neighborhoods and genomic features selected by _reservoir_featureSelection.R_

_reservoirPredict_PN.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the reservoir host of each virus using phylogenetic neighborhoods

_reservoirPredict_selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the reservoir host of each virus using genomic features selected by _reservoir_featureSelection.R_

_vectorPredict_featureSelection.R_ Uses gradient boosting machines in h2o to estimate average feature importances for predicting reservoir hosts across different training sets 

_vectorPredict_PN+selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the vector of each virus using phylogenetic neighborhoods and genomic features selected by _vectorPredict_featureSelection.R_

_vectorPredict_PN.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the vector of each virus using phylogenetic neighborhoods

_vectorPredict_selGen.R_ Uses gradient boosting machines in h2o to train and generate study wide (bagged) predictions of the vector of each virus using genomic features selected by _vectorPredict_featureSelection.R_

### Python script
_algo_comparison.py_ Compares the predictive power of a variety of competing machine learning algorithms to predict reservoir hosts, arthropod-borne transmission and vector taxa from all possible genomic features
