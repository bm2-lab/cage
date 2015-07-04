from __future__ import division
from pyfasta import Fasta
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency

def __DetectMh(lst_pmh):
    idx_csite = 11
    idx_idbeg = 14
    idx_idend = 15
    idx_mhseq = 21
    
    str_mhtype = ""
    int_mhlen = 0
    str_mhseq = 'None'
    int_csite = int(lst_pmh[idx_csite])
    int_idbeg = int(lst_pmh[idx_idbeg]) + 1
    int_idend = int(lst_pmh[idx_idend]) - 1

    str_pid = lst_pmh[idx_mhseq][0:len(lst_pmh[idx_mhseq])//2]
    str_id = lst_pmh[idx_mhseq][len(lst_pmh[idx_mhseq])//2:]
    for i in range(1, len(str_pid) + 1):
        if str_pid[-i] == str_id[-i]:
            int_mhlen += 1
        else:
            break
    if int_mhlen > 0:
        str_mhseq = str_id[-int_mhlen:]
        if not int_idbeg <= int_csite <= (int_idend - int_mhlen):
            int_mhlen = 0
    if int_mhlen < 5:
        str_mhtype = 'NHEJ'
    elif int_mhlen <= 30:
        str_mhtype = 'MMEJ'
    else:
        str_mhtype = 'SSA'
    lst_pmh.append(str_mhtype)
    lst_pmh.append(str(int_mhlen))
    lst_pmh.append(str_mhseq)
    lst_mh = lst_pmh
    return lst_mh



def ProcessMh(str_f_samind, str_refgem, str_of_mh):
    f_samind = open(str_f_samind, 'r')
    f_mh = open(str_of_mh, 'w')

    gn_samind = (str_line.strip().split('\t') for str_line in f_samind if str_line.strip() != '')
    idx_chr = 2
    idx_idtype = 13
    idx_idbeg = 14
    idx_idend = 15
    idx_idlen = 16
    gn_pmh = (lst_samind for lst_samind in gn_samind if lst_samind[idx_idtype] == 'sdel')

    str_refpath = '%s/%s.fa'% (os.environ['FASTADB'], str_refgem)
    fa_ref = Fasta(str_refpath)

    int_pmhbeg = 0
    int_pmhend = 0
    str_pmhseq = ''
    lst_mh = []
    for lst_pmh in gn_pmh:
        int_pmhbeg = int(lst_pmh[idx_idbeg]) + 1 - int(lst_pmh[idx_idlen])
        int_pmhend = int(lst_pmh[idx_idend]) - 1
        str_pmhseq = fa_ref[lst_pmh[idx_chr]][int_pmhbeg - 1:int_pmhend].upper()
        lst_pmh.append(str_pmhseq)
        lst_mh = __DetectMh(lst_pmh)
        f_mh.write('\t'.join(lst_mh) + '\n')
    f_mh.close()


def AnalyzeMh(str_f_mh, str_of_mnst):
    dfm_mh = pd.read_csv(str_f_mh, sep='\t', header=None)
    dfm_mh.columns = ['readID', 'batchID', 'chrom', 'sbeg', 'send', 'sgID', 'sgstrand', 'gene', 'sgbeg', 'sgend', 'sgseq', 'c_site', 'CIGAR', 'idtype', 'idbeg', 'idend', 'idlen', 'fm_status', 'factor', 'count', 'freq', 'premhseq', 'mhtype', 'mhlen', 'mhseq']

    func_mh = lambda x: pd.DataFrame(dict(batchID=x.batchID.unique(),
                                          sgID=x.sgID.unique(),
                                          allind=np.int(np.rint(x.freq.sum())),
                                          inf=np.int(np.rint(x.freq[x.fm_status=='INF'].sum())),
                                          otf=np.int(np.rint(x.freq[x.fm_status=='OTF'].sum())),
                                          mh=np.int(np.rint(x.freq[x.mhlen!=0].sum())),
                                          mh_inf=np.int(np.rint(x.freq[(x.mhlen!=0) & (x.fm_status=='INF')].sum())),
                                          mh_otf=np.int(np.rint(x.freq[(x.mhlen!=0) & (x.fm_status=='OTF')].sum())),
                                          NHEJ=np.int(np.rint(x.freq[x.mhtype=='NHEJ'].sum())),
                                          MMEJ=np.int(np.rint(x.freq[x.mhtype=='MMEJ'].sum())),
                                          r_mh=x.freq[x.mhlen!=0].sum()/x.freq.sum(),
                                          r_mmej=x.freq[x.mhtype=='MMEJ'].sum()/x.freq.sum()))

    dfm_mnst = dfm_mh.groupby(['batchID', 'sgID'], group_keys=False).apply(func_mh).reset_index(drop=True).reindex(columns=['batchID', 'sgID', 'allind', 'inf', 'otf', 'mh', 'mh_inf', 'mh_otf', 'NHEJ', 'MMEJ', 'r_mh', 'r_mmej'])

    func_mn = lambda x: pd.DataFrame(dict(sgID=x.sgID.unique(),
                                          allind=x.allind.sum(),
                                          inf=x.inf.sum(),
                                          otf=x.otf.sum(),
                                          mh=x.mh.sum(),
                                          mh_inf=x.mh_inf.sum(),
                                          nonmh_inf = x.inf.sum() - x.mh_inf.sum(),
                                          mh_otf=x.mh_otf.sum(),
                                          nonmh_otf = x.otf.sum() - x.mh_otf.sum(),
                                          NHEJ=x.NHEJ.sum(),
                                          MMEJ=x.MMEJ.sum(),
                                          r_mh=x.mh.sum()/x.allind.sum(),
                                          r_mmej=x.MMEJ.sum()/x.allind.sum()))
    dfm_mnst = dfm_mnst.groupby(['sgID'], group_keys=False).apply(func_mn).reset_index(drop=True).reindex(columns=['sgID', 'allind', 'inf', 'otf', 'mh', 'mh_inf', 'nonmh_inf', 'mh_otf', 'nonmh_otf', 'NHEJ', 'MMEJ', 'r_mh', 'r_mmej'])

    func_idp_test = lambda x: fisher_exact(np.array([[x.mh_otf, x.mh_inf], [x.nonmh_otf, x.nonmh_inf]]))[1] if np.any(np.array([[x.mh_otf, x.mh_inf], [x.nonmh_otf, x.nonmh_inf]]) < 5) else chi2_contingency(np.array([[x.mh_otf, x.mh_inf], [x.nonmh_otf, x.nonmh_inf]]))[1] 
    dfm_mnst['pv_mh_fs'] = dfm_mnst.apply(func_idp_test, axis=1)
    
    dfm_mnst.to_csv(str_of_mnst, sep='\t', index=None)
