# German Dataset: Independent Modeling

import pandas as pd
from hiclass import LocalClassifierPerNode
from copy import deepcopy
from sklearn.model_selection import train_test_split
from utils_ml import *

X = pd.read_csv("data/gcnlgg.csv", index_col=0)
X, Y = X.iloc[:, 5:], X["dummy"]
XTRAIN, XTEST, YTRAIN, YTEST = train_test_split(X, Y, test_size=PROP_TE, random_state=SEED)

# ----------------- Experimentation -----------------
lrc_holdout, lrc_cv = expt2x2(LRC, CV_CUSTOM, XTRAIN, YTRAIN, XTEST, YTEST)

gsvm_holdout, gsvm_cv = expt2x2(GSVM, CV_CUSTOM, XTRAIN, YTRAIN, XTEST, YTEST)

with pd.ExcelWriter("results_py/gcn_expts.xlsx") as ew:
    lrc_holdout.to_excel(ew, sheet_name="lrc_holdout")
    lrc_cv.to_excel(ew, sheet_name="lrc_cv")
    gsvm_holdout.to_excel(ew, sheet_name="gsvm_holdout")
    gsvm_cv.to_excel(ew, sheet_name="gsvm_cv")

# ----------------- Retrain Final Models & Export for Generalizability Analysis -----------------
lrc_flat = deepcopy(LRC)
lrc_flat.fit(X, Y)
save_model_as_pkl(lrc_flat, "./results_py/gcn_flat_lrc.pkl")

lrc_hc = LocalClassifierPerNode(local_classifier=deepcopy(LRC))
lrc_hc.fit(X, label_wrapper(Y))
save_model_as_pkl(lrc_hc, "./results_py/gcn_hc_lrc.pkl")

gsvm_flat = deepcopy(GSVM)
gsvm_flat.fit(X, Y)
save_model_as_pkl(gsvm_flat, "./results_py/gcn_flat_svm.pkl")

gsvm_hc = LocalClassifierPerNode(local_classifier=deepcopy(GSVM))
gsvm_hc.fit(X, label_wrapper(Y))
save_model_as_pkl(gsvm_hc, "./results_py/gcn_hc_svm.pkl")
