# Generalizability Analysis: American models on German cohort

import pandas as pd
from utils_ml import label_wrapper, apply_model_from_path

Xgcn = pd.read_csv("data/gcnlgg.csv", index_col=0)
Xgcn, Ygcn = Xgcn.iloc[:, 5:], Xgcn["dummy"]

res_gcn = {}
res_gcn["lrc_flat"] = apply_model_from_path("./results_py/tcga_flat_lrc.pkl", Xgcn, Ygcn)
res_gcn["lrc_hc"] = apply_model_from_path("./results_py/tcga_hc_lrc.pkl", Xgcn, label_wrapper(Ygcn), True)
res_gcn["svm_flat"] = apply_model_from_path("./results_py/tcga_flat_svm.pkl", Xgcn, Ygcn)
res_gcn["svm_hc"] = apply_model_from_path("./results_py/tcga_hc_svm.pkl", Xgcn, label_wrapper(Ygcn), True)
res_gcn = pd.DataFrame(res_gcn).transpose()
res_gcn.to_csv("./results_py/american2german.csv", index=True)