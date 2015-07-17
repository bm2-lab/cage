from __future__ import division
import numpy as np
import pandas as pd
from lxml import etree
from collections import namedtuple
from seqfeature_extractor import ExtractSeqFeature
from src.core import ml

Fr = namedtuple('Fr', ['fs', 'ups', 'dws'])

def __SaveFeatureReport(fr, str_of_fesrep):
    f_fesrep = open(str_of_fesrep, 'w')
    enode_root = etree.Element('report')
    enode_fes =etree.SubElement(enode_root, 'features')
    for fs in fr.fs:
        enode_fes.append(etree.Element(fs))
    enode_fes.set('method', 'mtlasso')
    enode_fes.set('n', str(len(fr.fs)))
    enode_fes.set('ups', str(fr.ups))
    enode_fes.set('dws', str(fr.dws))
    
    f_fesrep.write(etree.tostring(enode_root, pretty_print=True, xml_declaration=True, encoding='utf-8'))
    f_fesrep.close()

def MtLassoFeatureSelection(lst_sg, lst_st, lst_ref, lst_seq,
                            str_of_fesrep, int_ups=30, int_dws=27, flt_lmd=0.5):
    lst_x = []
    lst_y = []
    ind = [0]
    lst_header = []
    for i in xrange(len(lst_sg)):
        dfm = ExtractSeqFeature(lst_sg[i], lst_ref[i], int_ups, int_dws)
        dfm.to_csv(lst_seq[i], sep='\t', index=None)
        dfm_y = pd.read_csv(lst_st[i], sep='\t', index_col=None)[[0,-1]]
        dfm = pd.merge(dfm, dfm_y, on='sgID')
        x = np.array(dfm.ix[:,1:-1], dtype=np.double)
        y = np.array(dfm[[-1]], dtype=np.double)
        lst_x.append(x)
        lst_y.append(y)
        ind.append(dfm.shape[0])
        if i == 0:
            lst_header = dfm.columns[1:-1].tolist()
            
    x = np.double(np.vstack(lst_x))
    y = np.double(np.vstack(lst_y))
    ind = np.cumsum(ind)
    dict_opts = dict(init=2, rFlag=1, rsL2=0, ind=ind, nFlag=2)
    idx = ml.MtLassoSelector(x, y, flt_lmd, dict_opts)
    lst_fs = [lst_header[i] for i in idx]
    fr = Fr(fs=lst_fs, ups=int_ups, dws=int_dws)
    __SaveFeatureReport(fr, str_of_fesrep)
