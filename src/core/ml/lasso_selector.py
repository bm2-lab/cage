from __future__ import division
import numpy as np
from sklearn import linear_model
from sklearn.cross_validation import ShuffleSplit
from sklearn import preprocessing as prep
from sklearn import metrics
from sklearn.externals import joblib
from collections import namedtuple

Cvset = namedtuple('Cv', ['xtr', 'ytr', 'xte', 'yte'])
Md = namedtuple('Md', ['model', 'idx', 'r2'])

def LassoSelector(x, y, cv, niter):
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

    lr.fit(x[:,columns], y)
        
    return Md(model=lr.coef_, idx=columns, r2=np.mean(l1r2))

def SaveModel(model, str_of_md):
    joblib.dump(model, str_of_md)
