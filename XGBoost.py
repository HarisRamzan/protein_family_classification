# -*- coding: utf-8 -*-
"""
Created on Thu May  5 02:30:23 2020

@author: haris
"""

from warnings import simplefilter
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split,cross_val_score
import xgboost

def featurescaling(X_train,X_test):
    sc = MinMaxScaler(feature_range=[0, 1])
    #sc = StandardScaler()
    TrainX_std = sc.fit_transform(X_train)
    Testx_std = sc.fit_transform(X_test)
    return TrainX_std,Testx_std

def calculateXGBoostAccuracy(y_test,prediction):
    print("XGBoost accuracy:",accuracy_score(y_test,prediction))
    
def main():
    simplefilter(action='ignore', category=FutureWarning)
    #reading csv files
    Dataset= pd.read_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\FeatureExtracted.csv')
    Dataset =Dataset[Dataset['classification'] != 'Others'] #Dataset for binary classification
    print(Dataset)
    Dataset = Dataset.replace([np.inf, -np.inf,np.NaN],0)
    Y=Dataset['classification']
    X = Dataset.drop('classification', axis=1)
    X = X.drop('sequence', axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.3,random_state=5)
    print(X_train)
    TrainX_std,Testx_std=featurescaling(X_train,X_test);
    classifier = xgboost.XGBClassifier()
    classifier.fit(TrainX_std,y_train)
    prediction = classifier.predict(Testx_std)
    calculateXGBoostAccuracy(y_test, prediction)
    print("F score:",f1_score(y_test, prediction, average='macro'))
if __name__ == '__main__':
        main()