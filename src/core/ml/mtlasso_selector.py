import numpy as np
import slep


def MtLassoSelector(x, y, lmd, opts):
    coef, __, __ = slep.mtLeastR(x, y, lmd, opts)
    columns = np.where(coef[:,0]!=0)[0]
    return columns
