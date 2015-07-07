from __future__ import division
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.cross_validation import ShuffleSplit
from sklearn import preprocessing as prep
from sklearn import metrics
from collections import namedtuple
from lxml import etree

Cvset = namedtuple('Cv', ['xtr', 'ytr', 'xte', 'yte'])
Fitting = namedtuple('Fitting', ['under', 'over'])

def LassoFeatureSelection(str_f_seq, str_f_iost, str_of_fesrep, cv=10, niter=1000):
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

    l1r2 = Fitting(under=np.array([]), over=np.array([]))

    gn_cvsetl1 = (Cvset(func_prep(x[i][:,columns]), y[i], func_prep(x[j][:,columns]), y[j]) for (i, j) in ShuffleSplit(len(y), n_iter=niter, test_size=t_size))

    for cvt in gn_cvsetl1:
        lr.fit(cvt.xtr, cvt.ytr)
        l1r2 = l1r2._replace(under=np.hstack((l1r2.under, func_score(cvt.ytr, lr.predict(cvt.xtr)))), over=np.hstack((l1r2.over, func_score(cvt.yte, lr.predict(cvt.xte)))))
    lmn = (l1r2.under.mean(), l1r2.over.mean())

    int_ups = int(dfm_x.columns[1].split('_')[1])
    int_dws = int(dfm_x.columns[-1].split('_')[1])
    lst_fs = [i for i in dfm_x.columns[columns]]
    f_fesrep = open(str_of_fesrep, 'w')
    enode_root = etree.Element('report')
    enode_fes =etree.SubElement(enode_root, 'features')
    enode_cv = etree.SubElement(enode_root, 'cross_validation')
    for fs in lst_fs:
        enode_fes.append(etree.Element(fs))
    enode_fes.set('n', str(len(lst_fs)))
    enode_fes.set('ups', str(int_ups))
    enode_fes.set('dws', str(int_dws))
    enode_cv.set('n_iter', str(niter))
    enode_cv.set('fold', str(cv))
    enode_cv.set('r2', str(round(lmn[1], 3)))
    f_fesrep.write(etree.tostring(enode_root, pretty_print=True, xml_declaration=True, encoding='utf-8'))
    f_fesrep.close()
    

