import os
from indel_integrator import *
from src.core import seqfs

def AnalyzeIndel(opts):
    str_nm = os.path.basename(os.path.splitext(opts.samind)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_samind = opts.samind
    str_path_iost = os.path.join(str_path_proj, str_nm + '.iost')
    str_path_iosg = os.path.join(str_path_proj, str_nm + '.iosg')
    str_path_seq = os.path.join(str_path_proj, str_nm + '.seq')
    str_path_fesrep = os.path.join(opts.tdir, str_nm + '_fesrep.xml')
    str_proj = 'model'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_model = os.path.join(str_path_proj, str_nm + '.pkl')
    print('Integrating Indel information...')
    AnalyzeSamind(str_path_samind, opts.sg, str_path_iost, str_path_iosg, opts.cut)
    print('Done')
    if opts.auto == False:
        print('Selecting sequence feature...')
        seqfs.LassoFeatureSelection(str_path_iosg, str_path_iost, opts.ref, 
                              str_path_seq, str_path_fesrep, str_path_model,
                              opts.ups, opts.dws, opts.cv, opts.niter)
        print('Done')
    else:
        print('Detecting Optimal Feature Space...')
        seqfs.LassoRegionOptimizer(str_path_iosg, str_path_iost, opts.ref,
                             str_path_seq, str_path_fesrep, str_path_model,
                             opts.ir, opts.ir+opts.rad, opts.step,
                             opts.cv, opts.niter)
        print('Done')


