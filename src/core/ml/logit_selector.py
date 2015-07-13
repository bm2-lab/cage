from __future__ import division

import numpy as np
from sklearn import preprocessing as prep
from sklearn.cross_validation import ShuffleSplit
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

Cvset = namedtuple('Cv', ['xtr', 'ytr', 'xte', 'yte'])
Mdc = namedtuple('Mdc', ['model', 'idx', 'accu', 'prec', 'rec', 'f1', 'au'])

def __Auc(cls, xte, yte):
    ypo = cls.predict_proba(xte)
    flt_auc = roc_auc_score(yte, ypo[:,1])
    return flt_auc

def LogitSelector(x, y, cv, niter, njob):
    t_size=1 / cv

    lb = prep.LabelBinarizer()
    y = lb.fit_transform(y).ravel()

    model = LogisticRegressionCV(penalty='l1', solver='liblinear', refit=False, cv=cv, n_jobs=njob)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        warnings.simplefilter('ignore', ConvergenceWarning)
        model.fit(x, y)
    columns = np.arange(x.shape[1])[model.coef_.ravel() != 0]

    accu = []
    prec = []
    rec = []
    f1 = []
    au = []
    cls = LogisticRegression()
    gn_cvset = (Cvset(x[i][:, columns], y[i], x[j][:, columns], y[j]) for (i, j) in ShuffleSplit(len(y), n_iter=niter, test_size=t_size))

    for cvt in gn_cvset:
        cls.fit(cvt.xtr, cvt.ytr)
        accu.append(accuracy_score(cvt.yte, cls.predict(cvt.xte)))
        prec.append(precision_score(cvt.yte, cls.predict(cvt.xte)))
        rec.append(recall_score(cvt.yte, cls.predict(cvt.xte)))
        f1.append(f1_score(cvt.yte, cls.predict(cvt.xte)))
        au.append(__Auc(cls, cvt.xte, cvt.yte))

    cls.fit(x[:,columns], y)
    return Mdc(model=cls, idx=columns, accu=np.mean(accu),
               prec=np.mean(prec), rec=np.mean(rec), f1=np.mean(f1),
               au=np.mean(au))

