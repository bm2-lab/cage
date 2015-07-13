from __future__ import division
import numpy as np
import pandas as pd
from lxml import etree
from collections import namedtuple
from seqfeature_extractor import ExtractSeqFeature
from src.core import ml



Frc = namedtuple('Frc', ['fs', 'idx', 'ups', 'dws', 'niter', 'cv', 'accu', 'prec', 'rec', 'f1', 'au'])

def __SaveFeatureReport(frc, str_of_fesrep):
    f_fesrep = open(str_of_fesrep, 'w')
    enode_root = etree.Element('report')
    enode_fes =etree.SubElement(enode_root, 'features')
    enode_cv = etree.SubElement(enode_root, 'cross_validation')
    enode_met = etree.SubElement(enode_cv, 'metric')
    for (fs, idx) in zip(frc.fs, frc.idx):
        enode_fes.append(etree.Element(fs))
    enode_fes.set('method', 'logit')
    enode_fes.set('n', str(len(frc.fs)))
    enode_fes.set('ups', str(frc.ups))
    enode_fes.set('dws', str(frc.dws))
    
    enode_cv.set('n_iter', str(frc.niter))
    enode_cv.set('fold', str(frc.cv))

    enode_met.append(etree.Element('accuracy', value='%.3f'% frc.accu))
    enode_met.append(etree.Element('precision', value='%.3f'% frc.prec))
    enode_met.append(etree.Element('recall', value='%.3f'% frc.rec))
    enode_met.append(etree.Element('f1', value='%.3f'% frc.f1))
    enode_met.append(etree.Element('auc', value='%.3f'% frc.au))
    
    f_fesrep.write(etree.tostring(enode_root, pretty_print=True, xml_declaration=True, encoding='utf-8'))
    f_fesrep.close()


    
def LogitFeatureSelection(str_f_iosg, str_f_st, str_refgem,
                          str_of_seq, str_of_fesrep, str_of_md,
                          int_ups=30, int_dws=27, cv=5, niter=1000, njob=1):
    
    dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
    dfm_x.to_csv(str_of_seq, sep='\t', index=None)
    dfm_y = pd.read_csv(str_f_st, sep='\t', index_col=None)
    dfm_y = dfm_y[[0,-1]]
    dfm_y = pd.merge(dfm_x, dfm_y, on='sgID')
    dfm_x = dfm_y.ix[:, 1:-1]
    dfm_y = dfm_y[[-1]]

    x = np.array(dfm_x, dtype=np.double)
    y = np.array(dfm_y, dtype=np.int)

    mdc = ml.LogitSelector(x, y, cv, niter, njob)
    lst_fs = [i for i in dfm_x.columns[mdc.idx]]
    frc = Frc(fs=lst_fs, idx=mdc.idx, ups=int_ups, dws=int_dws, niter=niter,
              cv=cv, accu=mdc.accu, prec=mdc.prec, rec=mdc.rec,
              f1=mdc.f1, au=mdc.au)
    __SaveFeatureReport(frc, str_of_fesrep)
    mdlc = dict(med='logit', model=mdc.model, idx=mdc.idx,
                ups=int_ups, dws=int_dws)
    ml.SaveModel(mdlc, str_of_md)



def LogitRegionOptimizer(str_f_iosg, str_f_st, str_refgem,
                         str_of_seq, str_of_fesrep, str_of_md,
                         int_regn_start=0, int_regn_end=100, int_regn_step=5,
                         cv=5, niter=1000, njob=1):
    int_ups_regn = np.arange(int_regn_start, int_regn_end+1, int_regn_step)
    int_dws_regn = np.where((int_ups_regn-3)>0, int_ups_regn-3, 0)
    arr_accu = np.zeros_like(int_ups_regn, dtype=np.double)
    dfm_st = pd.read_csv(str_f_st, sep='\t', index_col=None)
    dfm_st = dfm_st[[0,-1]]
    for (i,(int_ups, int_dws)) in enumerate(zip(int_ups_regn, int_dws_regn)):
        dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
        dfm_y = pd.merge(dfm_x, dfm_st, on='sgID')
        dfm_x = dfm_y.ix[:, 1:-1]
        dfm_y = dfm_y[[-1]]
        x = np.array(dfm_x, dtype=np.double)
        y = np.array(dfm_y, dtype=np.int)
        arr_accu[i] = ml.LogitSelector(x, y, cv, niter, njob).accu
    int_ups = int_ups_regn[arr_accu.argmax()]
    int_dws = int_dws_regn[arr_accu.argmax()]
    dfm_x = ExtractSeqFeature(str_f_iosg, str_refgem, int_ups, int_dws)
    dfm_y = pd.merge(dfm_x, dfm_st, on='sgID')
    dfm_x = dfm_y.ix[:, 1:-1]
    dfm_y = dfm_y[[-1]]
    x = np.array(dfm_x, dtype=np.double)
    y = np.array(dfm_y, dtype=np.int)
    mdc = ml.LogitSelector(x, y, cv, niter, njob)
    lst_fs = [i for i in dfm_x.columns[mdc.idx]]
    frc = Frc(fs=lst_fs, idx=mdc.idx, ups=int_ups, dws=int_dws, niter=niter,
              cv=cv, accu=mdc.accu, prec=mdc.prec, rec=mdc.rec,
              f1=mdc.f1, au=mdc.au)
    __SaveFeatureReport(frc, str_of_fesrep)
    mdlc = dict(med='logit', model=mdc.model, idx=mdc.idx,
                ups=int_ups, dws=int_dws)
    ml.SaveModel(mdlc, str_of_md)
    

