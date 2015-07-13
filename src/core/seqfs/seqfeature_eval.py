import numpy as np
import pandas as pd
from lxml import etree
from collections import namedtuple
from seqfeature_extractor import ExtractSeqFeature
from src.core import ml


def FeatureEval(str_f_iosg, str_f_md, str_refgem, str_of_seq, str_of_st):
    mdl = ml.LoadModel(str_f_md)
    dfm = ExtractSeqFeature(str_f_iosg, str_refgem, mdl['ups'], mdl['dws'])
    dfm.to_csv(str_of_seq, sep='\t', index=None)

    x = np.array(dfm.ix[:,1:], dtype=np.double)
    y = None
    if mdl['med'] == 'lasso':
        y = ml.LassoEval(mdl['model'], x, mdl['idx'])
    elif mdl['med'] == 'logit':
        y = ml.LogitEval(mdl['model'], x, mdl['idx'])
    dfm['score'] = y

    dfm_y = dfm[[0,-1]]
    dfm_y.to_csv(str_of_st, sep='\t', index=None)
