# Utility Module for Flat ML & Hi-Class Analysis

import numpy as np
import pandas as pd
from copy import deepcopy
import pickle
from hiclass import LocalClassifierPerNode
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import cross_validate, KFold
from sklearn.metrics import make_scorer, confusion_matrix
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score

## Constants:
SEED = 240215
PROP_TE = 0.2

LRC = LogisticRegression(penalty="l1", solver="saga", multi_class="ovr")
GSVM = SVC(C=1, kernel="rbf", probability=True, decision_function_shape="ovr")

CV_CUSTOM = KFold(5, shuffle=True, random_state=SEED)
SCORERS_FLAT = ["accuracy","balanced_accuracy","precision_macro","recall_macro","f1_macro"]
SCORERS_HC = {
    "accuracy": make_scorer(lambda yt, yp: accuracy_score(undolabelwrapper(yt), undolabelwrapper(yp))),
    "balacc": make_scorer(lambda yt, yp: balanced_accuracy_score(undolabelwrapper(yt), undolabelwrapper(yp))),
    "precision": make_scorer(lambda yt, yp: precision_score(undolabelwrapper(yt), undolabelwrapper(yp), average="macro")),
    "recall": make_scorer(lambda yt, yp: recall_score(undolabelwrapper(yt), undolabelwrapper(yp), average="macro")),
    "f1_macro": make_scorer(lambda yt, yp: f1_score(undolabelwrapper(yt), undolabelwrapper(yp), average="macro"))
}

## Wrappers & helpers:
def expt2x2(learner, cvobj, Xtrain, ytrain, Xtest, ytest):
    """
    Wrapper to execute 2x2 experimental design
    :param learner: Un-fit sklearn multiclass model
    :return: 2 catcheable DataFrames
    """
    flat = deepcopy(learner)
    hcls = LocalClassifierPerNode(local_classifier=deepcopy(learner))

    ## Hold_out:
    df_holdout = {}

    flat.fit(Xtrain, ytrain)
    df_holdout["Flat"] = get_metrics(ytest, flat.predict(Xtest))

    hcls.fit(Xtrain, label_wrapper(ytrain))
    df_holdout["Hcls"] = get_metrics(ytest, undolabelwrapper(hcls.predict(Xtest)))

    df_holdout = pd.DataFrame(df_holdout)

    ## Kfold CV:
    Xall, yall = pd.concat([Xtrain, Xtest], axis=0), pd.concat([ytrain, ytest])
    cvscores_fl = cross_validate(flat, Xall, yall, cv=cvobj, scoring=SCORERS_FLAT)
    cvscores_hc = cross_validate(hcls, Xall, label_wrapper(yall), cv=cvobj, scoring=SCORERS_HC)

    df_cv = pd.DataFrame({
        "acc_flat": cvscores_fl["test_accuracy"],
        "acc_hcls": cvscores_hc["test_accuracy"],
        "balacc_flat": cvscores_fl["test_balanced_accuracy"],
        "balacc_hcls": cvscores_hc["test_balacc"],
        "precision_flat": cvscores_fl["test_precision_macro"],
        "precision_hcls": cvscores_hc["test_precision"],
        "recall_flat": cvscores_fl["test_recall_macro"],
        "recall_hcls": cvscores_hc["test_recall"],
        "F1_flat": cvscores_fl["test_f1_macro"],
        "F1_hcls": cvscores_hc["test_f1_macro"]
    })

    return df_holdout, df_cv

def apply_model_from_path(model_path, X, y, hcls=False):
    """
    Wrapper to apply models saved as pickle
    """
    model = pickle.load(open(model_path, "rb"))
    ypred = model.predict(X)
    if hcls:
        y = undolabelwrapper(y)
        ypred = undolabelwrapper(ypred)
    return get_metrics(y, ypred)

def label_wrapper(y):
    res = []
    for k in range(len(y)):
        if y[k] == 2:
            res.append(["WT"])
        elif y[k] == 1:
            res.append(["MUT", "C"])
        elif y[k] == 0:
            res.append(["MUT", "NC"])
    return np.array(res, dtype=object, ) #required syntax

def undolabelwrapper(y):
    res = []
    for k in range(len(y)):
        if y[k][0] == "WT":
            res.append(2)
        else:
            res.append(1 if y[k][1] == "C" else 0)
    return res

def get_metrics(yt, yp, average="macro"):
    print(confusion_matrix(yt, yp))
    return {
        "Accuracy": accuracy_score(yt, yp),
        "BalAcc": balanced_accuracy_score(yt, yp),
        "Precision": precision_score(yt, yp, average=average),
        "Recall": recall_score(yt, yp, average=average),
        "F1_"+average: f1_score(yt, yp, average=average)
    }

def save_model_as_pkl(clf, fname) -> None:
    with open(fname, 'wb') as f:
        pickle.dump(clf, f)
    print("Model saved at: %s" % fname)
