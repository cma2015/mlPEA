import sys
import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import train_test_split
from concurrent.futures import ThreadPoolExecutor

input_file = sys.argv[1]
output_file = sys.argv[2]
feature_num = int(sys.argv[3])
iteration = int(sys.argv[4])

df = pd.read_csv(input_file, sep="\t")

X = df.iloc[:, :-1]
y = df.iloc[:, -1]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

def compute_mi(feature):
    return mutual_info_classif(X_train[[feature]], y_train, discrete_features=True)[0]

with ThreadPoolExecutor() as executor:
    mi_scores = list(executor.map(compute_mi, X_train.columns))

mi_scores_df = pd.DataFrame(list(zip(X.columns, mi_scores)), columns=['Feature', 'Importance'])

mi_scores_df = mi_scores_df.sort_values(by='Importance', ascending=False)

top_features = mi_scores_df.head(feature_num)[['Feature', 'Importance']].values.tolist()

top_features_df = pd.DataFrame(top_features, columns=['Feature', 'Importance'])

output_file = f"{output_file}"
top_features_df.to_csv(output_file, sep=',', index=False, header=True)