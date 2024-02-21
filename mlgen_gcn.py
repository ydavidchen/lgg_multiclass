# Generalizability: German models on American cohort

import pandas as pd
from utils_ml import label_wrapper, apply_model_from_path

Xtcga = pd.read_csv("data/tcgalgg.csv", index_col=0)
Xtcga, Ytcga = Xtcga.iloc[:, 5:], Xtcga["dummy"]

res_tcga = {}
res_tcga["lrc_flat"] = apply_model_from_path("./results_py/gcn_flat_lrc.pkl", Xtcga, Ytcga)
res_tcga["lrc_hc"] = apply_model_from_path("./results_py/gcn_hc_lrc.pkl", Xtcga, label_wrapper(Ytcga), True)
res_tcga["svm_flat"] = apply_model_from_path("./results_py/gcn_flat_svm.pkl", Xtcga, Ytcga)
res_tcga["svm_hc"] = apply_model_from_path("./results_py/gcn_hc_svm.pkl", Xtcga, label_wrapper(Ytcga), True)
res_tcga = pd.DataFrame(res_tcga).transpose()
res_tcga.to_csv("./results_py/german2american.csv", index=True)