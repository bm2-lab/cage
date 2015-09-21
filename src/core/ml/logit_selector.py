from __future__ import division

import numpy as np
from sklearn import preprocessing as prep
from sklearn.cross_validation import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from collections import namedtuple
import warnings
from sklearn.utils import ConvergenceWarning

Mdc = namedtuple('Mdc', ['model', 'idx', 'accu', 'prec', 'rec', 'f1', 'au'])



def LogitSelector(x, y, cv, njob):

    lb = prep.LabelBinarizer()
    y = lb.fit_transform(y).ravel()

    cls = LogisticRegression()
    def __Auc(xte, yte):
        ypo = cls.predict_proba(xte)
        flt_auc = roc_auc_score(yte, ypo[:,1])
        return flt_auc
    
    skf = StratifiedKFold(y, n_folds=cv)
    model = LogisticRegressionCV(penalty='l1', solver='liblinear', fit_intercept=False, cv=cv, n_jobs=njob)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        warnings.simplefilter('ignore', ConvergenceWarning)
        model.fit(x, y)
    columns = np.arange(x.shape[1])[model.coef_.ravel() != 0]
    
    mdl_eval = lambda func: lambda idx_tr, idx_te: func(y[idx_te], cls.fit(x[idx_tr][:,columns], y[idx_tr]).predict(x[idx_te][:,columns]))
    auc_eval = lambda idx_tr, idx_te: roc_auc_score(y[idx_te], cls.fit(x[idx_tr][:,columns], y[idx_tr]).predict_proba(x[idx_te][:,columns])[:,1])
    res_eval = lambda func: np.average(map(mdl_eval(func), *zip(*[(idx_tr, idx_te) for idx_tr, idx_te in skf])))

    accu = res_eval(accuracy_score)
    prec = res_eval(precision_score)
    rec = res_eval(recall_score)
    f1 = res_eval(f1_score)
    au = np.average(map(auc_eval, *zip(*[(idx_tr, idx_te) for idx_tr, idx_te in skf])))

    cls.fit(x[:,columns], y)
    return Mdc(model=cls, idx=columns, accu=accu, prec=prec, rec=rec, f1=f1, au=au)

