'''
Babayan, Orton & Streicker

Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes

-- Algorithm comparison
'''

# COMPARISON BETWEEN ALGORITHMS USING DEFAULT SETTINGS

#%%
# modules and functions

import pandas as pd

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import cross_val_score

from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

import xgboost as xgb
from xgboost import XGBClassifier

import h2o
from h2o.estimators.gbm import H2OGradientBoostingEstimator

h2o.init()  # load H2O server
# h2o.cluster().shutdown() # uncomment if h2o clusters already running

# fix random seed for consistency between runs and algorithms
seed = 4

#############################################################################
# PREDICTING RESERVOIR TAXON
#############################################################################

# Load datasets from disk #

#%%
df_reservoirs = pd.read_csv("<filename of the appropriate database>")

df_h2o_reservoirs = h2o.import_file("<filename of the appropriate database>")

#%%
# cross-validation settings
validation_size = 0.3 # use 30% of the dataset for testing amd 70% for trainging
num_splits = 10 # repeat splitting 10 times
scoring = "accuracy" # a simple metric is sufficient
sss = StratifiedShuffleSplit(n_splits=num_splits, test_size=validation_size, random_state=seed) # take random splits preserving proportions of each class

# models for comparison
models = []
models.append(("KNN", KNeighborsClassifier()))
models.append(("LR", LogisticRegression(random_state=seed)))
models.append(("SGD", SGDClassifier(random_state=seed)))
models.append(("NB", GaussianNB()))
models.append(("SVM", SVC(random_state=seed)))
models.append(("RF", RandomForestClassifier(random_state=seed)))
models.append(("XGB", XGBClassifier(random_state=seed)))

h2o_spot = H2OGradientBoostingEstimator(
    nfolds=num_splits)

# label data
X = df_reservoirs
y = df_reservoirs.index

## for H2O models
X_h2o = df_h2o_reservoirs.columns
del X_h2o[0]
y_h2o = "class"

train, test = df_h2o_reservoirs.split_frame(ratios=[1-validation_size], seed=seed)

h2o_spot.train(x=X_h2o,
               y=y_h2o,
               training_frame=train,
               validation_frame=test)

h2o_cv_results = h2o_spot.cross_validation_metrics_summary(
).as_data_frame().iloc[0, 3:].values.astype("float64")

# prepare placeholders for results
results = []
names = []

# for Scikit models
for name, model in models:
    cv_results = cross_val_score(model, X, y, cv=sss, scoring=scoring)
    results.append(cv_results)
    names.append(name)
    msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
        name, cv_results.mean(), cv_results.std())
    print(msg)

# append H2O info
name = "H2O_GBM"
names.append(name)
results.append(h2o_cv_results)
msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
    name, h2o_cv_results.mean(), h2o_cv_results.std())
print(msg)


#############################################################################
# PREDICTING IF ARTHROPOD-BORNE
#############################################################################

# Load datasets from disk #
#%%
df_arth_borne = pd.read_csv("<filename of the appropriate database>")

df_h2o_arth_borne = h2o.import_file("<filename of the appropriate database>")


# Run model comparison #

#%%
# cross-val settings
validation_size = 0.3
num_splits = 5
scoring = "accuracy"
sss = StratifiedShuffleSplit(n_splits=num_splits, test_size=validation_size, random_state=seed)

# models for comparison
models = []
models.append(("KNN", KNeighborsClassifier()))
models.append(("LR", LogisticRegression()))
models.append(("SGD", SGDClassifier()))
models.append(("NB", GaussianNB()))
models.append(("SVM", SVC()))
models.append(("RF", RandomForestClassifier()))
models.append(("XGB", XGBClassifier()))

# data
X = df_arth_borne
y = df_arth_borne.index

## for H2O models
X_h2o = df_h2o_arth_borne.columns
y_h2o = "class"
X_h2o.remove(y_h2o)

df_h2o_arth_borne[0] = df_h2o_arth_borne[0].asfactor()
train, valid, test = df_h2o_arth_borne.split_frame(ratios=[0.7, 0.15], seed=seed)

h2o_spot = H2OGradientBoostingEstimator(nfolds=num_splits)

h2o_spot.train(x=X_h2o,
                y=y_h2o,
                training_frame=train,
                validation_frame=test)

h2o_cv_results = h2o_spot.cross_validation_metrics_summary().as_data_frame().iloc[0,3:].values.astype("float64")

# prepare placeholders for results
results = []
names = []

# for Scikit models
for name, model in models:
    cv_results = cross_val_score(model, X, y, cv=sss, scoring=scoring)
    results.append(cv_results)
    names.append(name)
    msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
        name, cv_results.mean(), cv_results.std())
    print(msg)

# append H2O info
name = "H2O_GBM"
names.append(name)
results.append(h2o_cv_results)
msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
    name, h2o_cv_results.mean(), h2o_cv_results.std())
print(msg)

#############################################################################
# PREDICTING VECTOR TAXON
#############################################################################

# Load datasets from disk #

#%%
df__vec_taxon = pd.read_csv("<filename of the appropriate database>")

df_h2o__vec_taxon = h2o.import_file("<filename of the appropriate database>")


# Run model comparison #

#%%
# cross-val settings
validation_size = 0.3
num_splits = 5
scoring = "accuracy"
sss = StratifiedShuffleSplit(n_splits=num_splits, test_size=validation_size, random_state=seed)

# models for comparison
models = []
models.append(("KNN", KNeighborsClassifier()))
models.append(("LR", LogisticRegression()))
models.append(("SGD", SGDClassifier()))
models.append(("NB", GaussianNB()))
models.append(("SVM", SVC()))
models.append(("RF", RandomForestClassifier()))
models.append(("XGB", XGBClassifier()))

h2o_spot = H2OGradientBoostingEstimator(
    nfolds=num_splits)

# data
X = df__vec_taxon
y = df__vec_taxon.index

## for H2O models
X_h2o = df_h2o__vec_taxon.columns
del X_h2o[0]
y_h2o = "class"

train, test = df_h2o__vec_taxon.split_frame(ratios=[1-validation_size], seed=seed)

h2o_spot.train(x=X_h2o,
                y=y_h2o,
                training_frame=train,
                validation_frame=test)

h2o_cv_results = h2o_spot.cross_validation_metrics_summary().as_data_frame().iloc[0,3:].values.astype("float64")

# prepare placeholders for results
results = []
names = []

# for Scikit models
for name, model in models:
    cv_results = cross_val_score(model, X, y, cv=sss, scoring=scoring)
    results.append(cv_results)
    names.append(name)
    msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
        name, cv_results.mean(), cv_results.std())
    print(msg)

# append H2O info
name = "H2O_GBM"
names.append(name)
results.append(h2o_cv_results)
msg = "Cross val score for {0}: {1:.2%} ± {2:.2%}".format(
    name, h2o_cv_results.mean(), h2o_cv_results.std())
print(msg)