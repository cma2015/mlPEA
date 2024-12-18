import os
import sys
import time
import docopt
import pickle
import logging
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.datasets import make_moons

#from baggingPU import BaggingClassifierPU
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
#from sklearn.naive_bayes import GaussianNB
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.linear_model import LogisticRegression
#from sklearn.neural_network import MLPClassifier
from tqdm import tqdm

PU_input = sys.argv[1] #'Ath_PU_input.txt'
PU_output = sys.argv[2]  #"Ath_PU_output.txt"

df_raw = pd.read_csv(PU_input, sep="\t")
print(df_raw.authentic.value_counts())
print('Has null values', df_raw.isnull().values.any())
df_tmp = df_raw.copy()
NON_LBL = [c for c in df_tmp.columns if c not in ['authentic', 'final']]
data_Praw = df_tmp[df_tmp['authentic']==1]
data_Uraw = df_tmp[df_tmp['authentic']==0]
data_P = data_Praw[NON_LBL].values
data_U = data_Uraw[NON_LBL].values
NP = data_P.shape[0]
NU = data_U.shape[0]

np.random.seed(0)
T = 500
K = NP
train_label = np.zeros(shape=(NP+K,))
train_label[:NP] = 1.0
n_oob = np.zeros(shape=(NU,))
f_oob = np.zeros(shape=(NU, 2))
for i in tqdm(range(T)):
    # Bootstrap resample
    bootstrap_sample = np.random.choice(np.arange(NU), replace=True, size=K)
    # Positive set + bootstrapped unlabeled set
    data_bootstrap = np.concatenate((data_P, data_U[bootstrap_sample, :]), axis=0)
    ## Train model
    model = SVC(C=1.0, kernel='rbf', gamma='auto', decision_function_shape='ovo', probability=True)
    
    model.fit(data_bootstrap, train_label)
    # Index for the out of the bag (oob) samples
    idx_oob = sorted(set(range(NU)) - set(np.unique(bootstrap_sample)))
    # Transductive learning of oob samples
    f_oob[idx_oob] += model.predict_proba(data_U[idx_oob])
    n_oob[idx_oob] += 1

predict_proba = f_oob[:, 1]/n_oob

with open(PU_output, "w") as wfile:
    for ii in range(len(data_Uraw.index.values)):
        wfile.write("{},{}\n".format(data_Uraw.index.values[ii], 
                                     predict_proba[ii]))
