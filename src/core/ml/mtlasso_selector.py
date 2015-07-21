import numpy as np
import slep


def MtLassoSelector(x, y, lmd, opts):
    slep.LoadModule('mtLeastR')
    coef, __, __ = slep.mtLeastR.mtLeastR(x, y, lmd, opts)
    columns = np.where(coef[:,0]!=0)[0]
    return columns
