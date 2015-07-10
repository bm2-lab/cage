import os
from src.core import indel

def GeneralLas(opts):
    str_nm = os.path.basename(os.path.splitext(opts.label)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.mkdir(str_path_proj)
    str_path_lab = opts.label
    str_path_sg = opts.sg
    str_path_seq = os.path.join(str_path_proj, str_nm + '.seq')
    str_path_fesrep = os.path.join(opts.tdir, str_nm + '_fesrep.xml')
    str_proj = 'model'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.mkdir(str_path_proj)
    str_path_model = os.path.join(str_path_proj, str_nm + '.pkl')
    if opts.auto == False:
        print('Selecting sequence feature...')
        indel.LassoFeatureSelection(str_path_sg, str_path_lab, opts.ref, 
                              str_path_seq, str_path_fesrep, str_path_model,
                              opts.ups, opts.dws, opts.cv, opts.niter)
        print('Done')
    else:
        print('Detecting Optimal Feature Space...')
        indel.LassoRegionOptimizer(str_path_sg, str_path_lab, opts.ref,
                             str_path_seq, str_path_fesrep, str_path_model,
                             opts.ir, opts.ir+opts.rad, opts.step,
                             opts.cv, opts.niter)
        print('Done')


