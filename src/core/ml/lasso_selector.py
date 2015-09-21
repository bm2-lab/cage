from __future__ import division

import numpy as np
from sklearn import linear_model
from sklearn.cross_validation import KFold
from sklearn import preprocessing as prep
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from collections import namedtuple
import warnings
from sklearn.utils import ConvergenceWarning

Md = namedtuple('Md', ['model', 'idx', 'cor', 'r2'])

def LassoSelector(x, y, cv, njob):
    cor_score = lambda x, y: pearsonr(x, y)[0]
    
    lr = linear_model.LinearRegression(n_jobs=njob)
    skf = KFold(len(y), n_folds=cv)
    model = linear_model.LassoLarsCV(fit_intercept=False, cv=cv, n_jobs=njob)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        warnings.simplefilter('ignore', ConvergenceWarning)
        model.fit(x, y)
    columns = np.arange(x.shape[1])[model.coef_ != 0]
    
    mdl_eval = lambda func: lambda idx_tr, idx_te: func(y[idx_te], lr.fit(x[idx_tr][:,columns], y[idx_tr]).predict(x[idx_te][:,columns]))
    res_eval = lambda func: np.average(map(mdl_eval(func), *zip(*[(idx_tr, idx_te) for idx_tr, idx_te in skf])))

    l1r2 = res_eval(r2_score)
    l1cor = res_eval(cor_score)

    lr.fit(x[:,columns], y)
        
    return Md(model=lr, idx=columns, cor=l1cor, r2=l1r2)
