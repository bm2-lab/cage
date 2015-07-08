
from __future__ import division
import numpy as np
import pandas as pd

def AnalyzeSamind(str_f_samind, str_f_sg, str_of_iost, str_of_iosg, int_cutoff=0):
    dfm_samind = pd.read_csv(str_f_samind, sep='\t', header=None)
    dfm_samind.columns = ['readID', 'batchID', 'chrom', 'sbeg', 'send', 'sgID', 'sgstrand', 'gene', 'sgbeg', 'sgend', 'sgseq', 'c_site', 'CIGAR', 'idtype', 'idbeg', 'idend', 'idlen', 'fm_status', 'factor', 'count', 'freq']

    dfm_sg = pd.read_csv(str_f_sg, header=None, index_col=None, sep='\t')
    dfm_sg.columns = ['sgID', 'chrom', 'strand', 'gene', 'sbeg', 'send', 'qseq', 'c_site']

    func_samind = lambda x: pd.DataFrame(dict(batchID=x.batchID.unique(),
                                           sgID=x.sgID.unique(),
                                           allrc = np.int(np.rint(x.freq.sum())),
                                           none=np.int(np.rint(x.freq[x.fm_status=='None'].sum())),
                                              r_ind = x.freq[x.fm_status!='None'].sum()/x.freq.sum() if x.freq.sum()!=0 else 0,
                                           inf=np.int(np.rint(x.freq[x.fm_status=='INF'].sum())),
                                           otf=np.int(np.rint(x.freq[x.fm_status=='OTF'].sum())),
                                           r_reffect=x.freq[x.fm_status=='OTF'].sum()/x.freq[x.fm_status!='None'].sum() if x.freq[x.fm_status!='None'].sum()!=0 else 0,
                                              r_effect=x.freq[x.fm_status=='OTF'].sum()/x.freq.sum() if x.freq.sum()!=0 else 0))
    
    dfm_iost = dfm_samind.groupby(['batchID', 'sgID'], group_keys=False).apply(func_samind).reset_index(drop=True).reindex(columns=['batchID', 'sgID', 'allrc', 'none', 'r_ind', 'inf', 'otf', 'r_reffect', 'r_effect'])

    func_io = lambda x: pd.DataFrame(dict(sgID=x.sgID.unique(),
                                          allrc=x.allrc.sum(),
                                          none=x.none.sum(),
                                          inf=x.inf.sum(),
                                          otf=x.otf.sum(),
                                          r_ind=(x.inf.sum()+x.otf.sum())/x.allrc.sum() if x.allrc.sum()!=0 else 0,
                                          r_reffect=x.otf.sum()/(x.inf.sum()+x.otf.sum()) if x.inf.sum()+x.otf.sum()!=0 else 0,
                                          r_effect=x.otf.sum()/x.allrc.sum() if x.allrc.sum()!=0 else 0))
    dfm_iost = dfm_iost.groupby(['sgID'], group_keys=False).apply(func_io).reset_index(drop=True).reindex(columns=['sgID', 'allrc', 'none', 'inf', 'otf', 'r_ind', 'r_reffect', 'r_effect'])
    dfm_iost = dfm_iost.ix[dfm_iost.allrc>=int_cutoff, ]
    dfm_iost = dfm_iost.sort_index(by=['r_effect'], ascending=[False])
    dfm_iost.to_csv(str_of_iost, sep='\t', index=None)

    dfm_sg = pd.merge(dfm_sg, dfm_iost[['sgID']], on='sgID')
    dfm_sg.to_csv(str_of_iosg, sep='\t', header=None, index=None)

    
