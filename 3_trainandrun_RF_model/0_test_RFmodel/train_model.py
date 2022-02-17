#! /usr/local/bin/python3

# os.system("pip install pandas numpy scikit-learn scipy loompy p_tqdm")

import os
import yaml
import pickle
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score

import logging
logging.basicConfig()
from pprint import pprint, pformat


def main():
    
    # CONSTS
    loomfile = os.getenv("LOOMFILE")
    metafile = os.getenv("METAFILE")
    paramfile = os.getenv("PARAMFILE")
    outmodel = os.getenv("OUTMODEL")
    outresult_xtra = os.getenv("OUTRESULT_5FOLDCV")
    outresult = os.getenv("OUTRESULT_UPSAMPLE")
    classcolumn = os.getenv("CLASSCOLUMN")
    
    # READ DATA
    adata = sc.read_loom(loomfile,
                         X_name="counts",
                         sparse=True, validate=False)
    # read meta
    meta_data = pd.read_csv(metafile, sep="\t")
    meta_data.index = meta_data.CellID
    
    # read params
    with open(paramfile, "r") as fl:
        params = yaml.load(fl, Loader=yaml.FullLoader)

      # merge meta
    mrgdf = meta_data.loc[adata.obs.index, :]
    adata.obs = adata.obs.combine_first(mrgdf)
    
    # remove doublets
#     adata = adata[(adata.obs.is_doublet == False) & (adata.obs[classcolumn].isna() == False), :]
    adata = adata[adata.obs.type.str.contains("(Doub)|(LowQual)") == False, :]

     # fit and save trained model
    fit_rf_model_and_save(
        model=RandomForestClassifier(n_jobs=-1),
        data=adata.layers["norm_data"],
        metadata=adata.obs,
#         rows=adata.obs.index,
        classname=classcolumn,
        params=params,
        modelfile=outmodel,
        resultsfile=outresult_xtra, 
        traintestsplit=5,
        logginglevel="INFO"
    )

    # overwrite model with better upsampled train test split
    # fit and save trained model
    fit_rf_model_and_save(
        model=RandomForestClassifier(n_jobs=-1),
        data=adata.layers["norm_data"],
        metadata=adata.obs,
#         rows=adata.obs.index,
        classname=classcolumn,
        params=params,
        modelfile=outmodel,
        resultsfile=outresult, 
        traintestsplit=ensure_even_class_representation_traintest_split(
            X=adata.obs.reset_index(), class_col=classcolumn, nsampling=5,
            sample_n_from_each_class = int(np.quantile(adata.obs.loc[:, classcolumn].value_counts(), .95))),
        logginglevel="INFO"
    )




            
def ensure_even_class_representation_traintest_split(
    X, class_col, nsampling=5, initialtrainfractionperclass=0.8, sample_n_from_each_class = 200, randomseed=None):
    """ Splits dataset into train & test sets
    
    test set ensures `1-initialtrainfractionperclass=0.2` of each class is in test set
    train set ensures `initialtrainfractionperclass=0.8` of each class is in training set 
        and upsamples to `sample_n_from_each_class=200` samples per class to help remove uneven class bias
        
    params:
        X = data frame to select rows from
        class_col = name of class column in X
        nsampling = number of time to perform train test split
        initialtrainfractionperclass = fraction of each class to sample without replacement to put in training set 
            (rest are placed in test set)
        sample_n_from_each_class = number to sample with replacement to for all classes in training set
        randomseed = random seed to replicate splits
        
    returns:
        list (length nsampling long) of tuples: _[0] = training set index, _[1] = test set index
    
    """
    np.random.seed(randomseed)   
    N = X.shape[0]
    for f in range(nsampling):
        N = X.shape[0]
        # split to ensure representation of each class
        trnidx = X.groupby(class_col)\
            .apply(lambda x: x.sample(frac=initialtrainfractionperclass, replace=False))\
            .index.get_level_values(1)
        tstidx = X.drop(trnidx).index
        
        assert N == X.loc[trnidx].shape[0] + X.drop(trnidx).shape[0]
        
        trnidx = X.loc[trnidx]\
            .groupby(class_col)\
            .apply(lambda x: x.sample(n=sample_n_from_each_class, replace=True))\
            .index.get_level_values(1)
             
        yield trnidx, tstidx
        
        
def fit_rf_model_and_save(
    model, # sklearn model able to be put in GridSearchCV model
    data, # scipy.csr_matrix
    metadata, # pd.dataframe that matches data on row indices
    classname, # column in metadata to use as class names
    params, # params dict for grid search
    modelfile, # where to save model pickle file
    resultsfile, # where to save model results (CV: accuracy, etc.)
    traintestsplit=5, # if 5 performs basic 5-fold CV can provide generator of index labels for data (must be int type)
    logginglevel="DEBUG" # login level if you want quiet output set "ERROR"
    ): 
    """
    """
    
    fn_logger = logging.getLogger("fit_rf_model_and_save")
    fn_logger.setLevel(logginglevel)
    
    
    # optimize parameters with grid search
    fn_logger.info(f"Starting grid search optimization over:\n{pformat(params)}")
    
    scoring = {'accuracy': make_scorer(accuracy_score),
                'precision_macro': make_scorer(precision_score, average = 'macro'),
                'precision_micro': make_scorer(precision_score, average = 'micro'),
                'precision_weighted': make_scorer(precision_score, average = 'weighted'),
                # 'precision_all': make_scorer(precision_score, average = None),
                'recall_macro': make_scorer(recall_score, average = 'macro'),
                'recall_micro': make_scorer(recall_score, average = 'micro'),
                'recall_weighted': make_scorer(recall_score, average = 'weighted'),
                # 'recall_all': make_scorer(recall_score, average = None),
                'f1_macro': make_scorer(f1_score, average = 'macro'),
                'f1_micro': make_scorer(f1_score, average = 'micro'),
                'f1_weighted': make_scorer(f1_score, average = 'weighted')}

    # https://scikit-learn.org/stable/modules/model_evaluation.html#multilabel-ranking-metrics
    optimized = GridSearchCV(estimator=model,
                   param_grid=params,
                   refit="accuracy",
                   scoring=scoring,
                   cv=traintestsplit,
                   n_jobs=-1, verbose=10)
    
#     fn_logger.debug(f"data:\n{data}")
#     fn_logger.debug(f"metadata:\n{metadata.loc[:, classname].reset_index(drop=True)}")
    fn_logger.debug(f"optimized:\n{pformat(optimized)}")
    
    optimized = optimized.fit(
        data, metadata.loc[:, classname].reset_index(drop=True).fillna("NAN"))
    fn_logger.info("Finished grid search optimization")
    
    ### write out results ###
    results = pd\
        .DataFrame(optimized.cv_results_)\
        .sort_values("rank_test_accuracy")\
        .T
    results.to_csv(resultsfile)
    fn_logger.info(f"wrote results to: {resultsfile}")
    
    
    
    ### refit model on best params ### 
    
    fn_logger.info(f"Refitting for final model")
    final_model = RandomForestClassifier(n_jobs=-1, **optimized.best_params_)
    
    trnidx, _ = [x for x in ensure_even_class_representation_traintest_split(
        metadata.loc[:].reset_index(), classname,
        initialtrainfractionperclass = 1.0,
        sample_n_from_each_class = max(metadata.loc[:, classname].value_counts()), 
        nsampling=1
    )][0]
    
    final_model = final_model.fit(
        data[trnidx, :], metadata.loc[:, classname].reset_index(drop=True)[trnidx].fillna("NAN"))
    fn_logger.info(f"dumping file to {modelfile}")
    with open(modelfile, "wb") as mf:
        pickle.dump(final_model, mf)
    fn_logger.info(f"Completed job!")
    
    
if __name__ == '__main__':
    main()