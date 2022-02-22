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

from sklearn import datasets
import numpy as np
iris = datasets.load_iris()
iris.target = np.random.randint(0,2,len(iris.target))

X = iris.data
y = iris.target


factory = {'svm' : (SVC(kernel = 'rbf', 
                        gamma = 'scale', 
                        probability = True, 
                        cache_size = 500, 
                        class_weight = 'balanced'), 
                    
                    {'C' : [10**i for i in range(-3,5)]}),
           
           'for' : (RandomForestClassifier(
                                           class_weight = 'balanced',
                                           n_estimators = 5000),
               
                    {'min_samples_leaf' : [1,2,3,4,6,8,10,15]}),
           
           'log' : (LogisticRegression(solver = 'liblinear', 
                                       penalty = 'l2', 
                                       max_iter = 200, 
                                       class_weight = 'balanced'),
                    
                    {'C' : [10**i for i in range(-7,1)]})}

CV_ROUND = 5
result = {}
for nm in list(factory.keys()):
    print(nm)
    model = factory[nm][0]
    parameters = factory[nm][1]
    clf = GridSearchCV(model, parameters, cv=CV_ROUND)
    scores = cross_val_score(clf, X, y, cv=CV_ROUND, scoring='roc_auc')
    result[nm] = scores
    



    


