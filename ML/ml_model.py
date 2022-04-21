#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:42:53 2022

@author: jefft
"""
### Machine learning models with optimization

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score
import numpy as np
import joblib
from sklearn.metrics import f1_score, make_scorer, roc_auc_score
from sklearn.preprocessing import label_binarize
from sklearn.utils import shuffle


'''
from sklearn import datasets
iris = datasets.load_iris()
iris.target = np.random.randint(0,2,len(iris.target))
X = iris.data
y = iris.target
'''

'''
from pyreadr import pyreadr
X = pyreadr.read_r('X_top2k.rds')[None]
y = pyreadr.read_r('y.rds')[None]
y = np.array(list(y[None]))
X = np.array(X)
'''

def run_opt(X,y,model=None,cv_round=5,sample_weight=None):
    factory = {'svm' : (SVC(gamma = 'auto', 
                            kernel = 'linear',
                            probability = True, 
                            cache_size = 500, 
                            class_weight = 'balanced'), 
                        
                        {'C' : [10**i for i in range(-3,5)]}),
            
            'for' : (RandomForestClassifier(oob_score = False,
                                            class_weight = 'balanced', random_state=3180110984),
                
                        {'min_samples_leaf' : [4,8,16,32],
                         'max_depth': np.arange(6,12,2),
                         'n_estimators': np.arange(100,900,300)}),
            
            'log' : (LogisticRegression(solver = 'liblinear', 
                                        penalty = 'l2', 
                                        max_iter = 200, 
                                        class_weight = 'balanced'),
                        
                        {'C' : [10**i for i in range(-7,1)]})}

    CV_ROUND = int(cv_round)
    # result = {}
    del factory['log']
    # del factory['svm']
    del factory['for']
    y = np.array(y)
    X, y = shuffle(X, y, random_state=3180110984)
    # y = label_binarize(y, classes=np.unique(y))
    f1_scorer = make_scorer(f1_score, average='micro')
    roc_scorer = make_scorer(roc_auc_score, average='macro', multi_class='ovo')
    for nm in list(factory.keys()):
        # print(nm)
        model = factory[nm][0]
        parameters = factory[nm][1]
        clf = GridSearchCV(model, parameters, cv=CV_ROUND, n_jobs=7, verbose=4,
                          scoring=f1_scorer
                           )
        print(X.dtype,y.dtype, X.shape, y.shape, X.__class__, y.__class__)
        clf.fit(X, y, sample_weight=sample_weight)
        bp = clf.best_params_
        #scores = cross_val_score(clf, X, y, cv=CV_ROUND, scoring='balanced_accuracy')
        #result[nm] = scores
    return bp


def eval_cls(X,y,model_fp,cv_round=5,scoring=make_scorer(roc_auc_score, avarage='macro', multi_class='ovo')):
    clf = joblib.load(model_fp)
    X, y = shuffle(X,y, random_state=3180110984)
    y = label_binarize(y, classes=np.unique(y))
    scores = cross_val_score(clf, X, y, cv=cv_round, scoring=scoring)
    return scores


def run_fit_forest(X,y,params,save_fp=None,sample_weight=None):
    clf = RandomForestClassifier(oob_score = True, class_weight = 'balanced', 
                                 n_estimators = int(params['n_estimators']), 
                                 min_samples_leaf = int(params['min_samples_leaf']),
                                 max_depth = int(params['max_depth']), random_state=318011984)
    clf.fit(X, y, sample_weight=sample_weight)
    if save_fp:
        joblib.dump(clf, save_fp)
    return {'predict_mtx':clf.predict_proba(X), 'feature_importance':clf.feature_importances_, 'oob':clf.oob_score_}


def run_fit_SVM(X,y,params,save_fp=None, sample_weight=None):
    clf = SVC(class_weight = 'balanced', kernel = 'linear', gamma = 'auto', 
             probability = True, cache_size = 200, 
             C = float(params['C']), random_state=318011984)
    clf.fit(X, y, sample_weight=sample_weight)
    if save_fp:
        joblib.dump(clf, save_fp)
    return {'predict_mtx':clf.predict_proba(X), 'coef':clf.coef_}


def run_cls(X, model_fp):
    clf = joblib.load(model_fp)
    return clf.predict_proba(X)
    


