# -*- coding: utf-8 -*-
"""
Created on Thu Apr  24 00:18:07 2020

@author: haris
"""
#Implementing svm for sequence data
from warnings import simplefilter
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
from sklearn.svm import LinearSVC,SVC
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split,cross_val_score

def featurescaling(X_train,X_test):
    sc = MinMaxScaler(feature_range=[0, 1])
    #sc = StandardScaler()
    TrainX_std = sc.fit_transform(X_train)
    Testx_std = sc.fit_transform(X_test)
    return TrainX_std,Testx_std

def TrainLiniersvmModel(svm,TrainX_std,y_train):
        svm.fit(TrainX_std,y_train)
        
def Predictliniersvm(svm,TestX_std):
    y_pred=svm.predict(TestX_std)
    return y_pred

def calculateliniersvmAccuracy(y_test,prediction):
    print("svc accuracy:",accuracy_score(y_test,prediction))


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
    X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.2,random_state=10)
    print(X_train)
    TrainX_std,Testx_std=featurescaling(X_train,X_test);
    svm=LinearSVC(random_state=0, tol=1e-5)
#    svm =SVC(C=5, cache_size=200, class_weight=None, coef0=0.0,
#    decision_function_shape='ovr', degree=3, gamma=0.01, kernel='rbf',
#    max_iter=-1, probability=False, random_state=None, shrinking=True,
#    tol=0.001, verbose=False)
    TrainLiniersvmModel(svm,TrainX_std,y_train);
    prediction=Predictliniersvm(svm,Testx_std)
    calculateliniersvmAccuracy(y_test, prediction)
    print(f1_score(y_test, prediction, average='macro'))
if __name__ == '__main__':
        main()