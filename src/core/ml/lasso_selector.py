from __future__ import division

import numpy as np
from sklearn import linear_model
from sklearn.cross_validation import ShuffleSplit
from sklearn import preprocessing as prep
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from collections import namedtuple
import warnings
from sklearn.utils import ConvergenceWarning

Cvset = namedtuple('Cv', ['xtr', 'ytr', 'xte', 'yte'])
Md = namedtuple('Md', ['model', 'idx', 'cor', 'r2'])

def LassoSelector(x, y, cv, niter, njob):
    t_size=1 / cv
    cor_score = lambda x, y: pearsonr(x, y)[0]
    
    lr = linear_model.LinearRegression(n_jobs=njob)
    model = linear_model.LassoLarsCV(fit_intercept=False, cv=cv, n_jobs=njob)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        warnings.simplefilter('ignore', ConvergenceWarning)
        model.fit(x, y)
    columns = np.arange(x.shape[1])[model.coef_ != 0]

    l1cor = []
    l1r2 = []

    gn_cvsetl1 = (Cvset(x[i][:,columns], y[i], x[j][:,columns], y[j]) for (i, j) in ShuffleSplit(len(y), n_iter=niter, test_size=t_size))

    for cvt in gn_cvsetl1:
        lr.fit(cvt.xtr, cvt.ytr)
        l1cor.append(cor_score(cvt.yte, lr.predict(cvt.xte)))
        l1r2.append(r2_score(cvt.yte, lr.predict(cvt.xte)))

    lr.fit(x[:,columns], y)
        
    return Md(model=lr, idx=columns, cor=np.mean(l1cor), r2=np.mean(l1r2))
