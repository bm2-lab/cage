from __future__ import division
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.cross_validation import ShuffleSplit
from sklearn import preprocessing as prep
from sklearn import metrics
from collections import namedtuple
from seqfeature_extractor import ExtractSeqFeature
from seqfeature_selector import LassoFeatureSelection

Cvset = namedtuple('Cv', ['xtr', 'ytr', 'xte', 'yte'])

def __LassoFs(str_f_seq, str_f_iost, cv=10, niter=1000):
    dfm_s_ori = pd.read_csv(str_f_seq, sep='\t', index_col=None)
    dfm_ef_ori = pd.read_csv(str_f_iost, sep='\t', index_col=None)
    dfm_ef = dfm_ef_ori[['sgID', 'r_effect']]
    dfm = pd.merge(dfm_ef, dfm_s_ori, on='sgID')

    dfm_x = dfm.ix[:, 2:]
    dfm_y = dfm.r_effect

    x = np.array(dfm_x, dtype=np.double)
    y = np.array(dfm_y, dtype=np.double)

    t_size=1 / cv
    func_score = metrics.r2_score
    func_prep = lambda x: x

    lr = linear_model.LinearRegression(normalize=False, fit_intercept=True)
    model = linear_model.LassoLarsCV(fit_intercept=False, cv=cv)

    model.fit(x, y)
    columns = np.arange(x.shape[1])[model.coef_ != 0]

    l1r2 = []

    gn_cvsetl1 = (Cvset(func_prep(x[i][:,columns]), y[i], func_prep(x[j][:,columns]), y[j]) for (i, j) in ShuffleSplit(len(y), n_iter=niter, test_size=t_size))

    for cvt in gn_cvsetl1:
        lr.fit(cvt.xtr, cvt.ytr)
        l1r2.append(func_score(cvt.yte, lr.predict(cvt.xte)))

    return np.mean(l1r2)

def __DetectRegion(str_f_iosg, str_f_iost, str_refgem, str_of_seq, int_regn_start=0, int_regn_end=80, int_regn_step=5, cv=10, niter=1000):

    int_ups_regn = np.arange(int_regn_start, int_regn_end+1, int_regn_step)
    int_dws_regn = np.where((int_ups_regn-3)>0, int_ups_regn-3, 0)
    arr_r2 = np.zeros_like(int_ups_regn, dtype=np.double)
    for (i,(int_ups, int_dws)) in enumerate(zip(int_ups_regn, int_dws_regn)):
        ExtractSeqFeature(str_f_iosg, str_refgem, str_of_seq, int_ups, int_dws)
        arr_r2[i] = __LassoFs(str_of_seq, str_f_iost, cv, niter)
    return (int_ups_regn, int_dws_regn, arr_r2)

def FeatureSelection(str_f_iosg, str_f_iost, str_refgem, str_of_seq, str_of_fesrep, int_regn_start=0, int_regn_end=80, int_regn_step=5, cv=10, niter=1000):
    int_ups_regn, int_dws_regn, arr_r2 = __DetectRegion(str_f_iosg, str_f_iost, str_refgem, str_of_seq, int_regn_start, int_regn_end, int_regn_step, cv, niter)

    int_ups = int_ups_regn[arr_r2.argmax()]
    int_dws = int_dws_regn[arr_r2.argmax()]
    ExtractSeqFeature(str_f_iosg, str_refgem, str_of_seq, int_ups, int_dws)
    LassoFeatureSelection(str_of_seq, str_f_iost, str_of_fesrep, cv, niter)
