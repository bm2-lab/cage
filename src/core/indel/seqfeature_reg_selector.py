from __future__ import division
import numpy as np
import pandas as pd
from lxml import etree
from collections import namedtuple
from seqfeature_extractor import ExtractSeqFeature
from src.core import ml

Fr = namedtuple('Fr', ['fs', 'ups', 'dws', 'niter', 'cv', 'r2'])

def __SaveFeatureReport(fr, str_of_fesrep):
    f_fesrep = open(str_of_fesrep, 'w')
    enode_root = etree.Element('report')
    enode_fes =etree.SubElement(enode_root, 'features')
    enode_cv = etree.SubElement(enode_root, 'cross_validation')
    for fs in fr.fs:
        enode_fes.append(etree.Element(fs))
    enode_fes.set('n', str(len(fr.fs)))
    enode_fes.set('ups', str(fr.ups))
    enode_fes.set('dws', str(fr.dws))
    enode_cv.set('n_iter', str(fr.niter))
    enode_cv.set('fold', str(fr.cv))
    enode_cv.set('r2', str(round(fr.r2, 3)))
    f_fesrep.write(etree.tostring(enode_root, pretty_print=True, xml_declaration=True, encoding='utf-8'))
    f_fesrep.close()


    
def LassoFeatureSelection(str_f_iosg, str_f_iost, str_refgem,
                          str_of_seq, str_of_fesrep, str_of_md,
                          int_ups=30, int_dws=27, cv=5, niter=1000):
    
    dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
    dfm_x.to_csv(str_of_seq, sep='\t', index=None)
    dfm_y = pd.read_csv(str_f_iost, sep='\t', index_col=None)
    dfm_y = dfm_y[[0,-1]]
    dfm_y = pd.merge(dfm_x, dfm_y, on='sgID')
    dfm_x = dfm_y.ix[:, 1:-1]
    dfm_y = dfm_y[[-1]]

    x = np.array(dfm_x, dtype=np.double)
    y = np.array(dfm_y, dtype=np.double)

    md = ml.LassoSelector(x, y, cv, niter)
    lst_fs = [i for i in dfm_x.columns[md.idx]]
    fr = Fr(fs=lst_fs, ups=int_ups, dws=int_dws, niter=niter, cv=cv, r2=md.r2)
    __SaveFeatureReport(fr, str_of_fesrep)
    ml.SaveModel(md.model, str_of_md)



def LassoRegionOptimizer(str_f_iosg, str_f_iost, str_refgem,
                         str_of_seq, str_of_fesrep, str_of_md,
                         int_regn_start=0, int_regn_end=100, int_regn_step=5,
                         cv=5, niter=1000):
    int_ups_regn = np.arange(int_regn_start, int_regn_end+1, int_regn_step)
    int_dws_regn = np.where((int_ups_regn-3)>0, int_ups_regn-3, 0)
    arr_r2 = np.zeros_like(int_ups_regn, dtype=np.double)
    dfm_iost = pd.read_csv(str_f_iost, sep='\t', index_col=None)
    dfm_iost = dfm_iost[[0,-1]]
    for (i,(int_ups, int_dws)) in enumerate(zip(int_ups_regn, int_dws_regn)):
        dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
        dfm_y = pd.merge(dfm_x, dfm_iost, on='sgID')
        dfm_x = dfm_y.ix[:, 1:-1]
        dfm_y = dfm_y[[-1]]
        x = np.array(dfm_x, dtype=np.double)
        y = np.array(dfm_y, dtype=np.double)
        arr_r2[i] = ml.LassoSelector(x, y, cv, niter).r2
    int_ups = int_ups_regn[arr_r2.argmax()]
    int_dws = int_dws_regn[arr_r2.argmax()]
    dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
    dfm_y = pd.merge(dfm_x, dfm_iost, on='sgID')
    dfm_x = dfm_y.ix[:, 1:-1]
    dfm_y = dfm_y[[-1]]
    x = np.array(dfm_x, dtype=np.double)
    y = np.array(dfm_y, dtype=np.double)
    md = ml.LassoSelector(x, y, cv, niter)
    lst_fs = [i for i in dfm_x.columns[md.idx]]
    fr = Fr(fs=lst_fs, ups=int_ups, dws=int_dws, niter=niter, cv=cv, r2=md.r2)
    __SaveFeatureReport(fr, str_of_fesrep)
    ml.SaveModel(md.model, str_of_md)
    

