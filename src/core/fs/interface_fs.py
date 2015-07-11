import os
from src.core import seqfs

def GeneralFs(opts):
    str_nm = os.path.basename(os.path.splitext(opts.label)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_lab = opts.label
    str_path_sg = opts.sg
    str_path_seq = os.path.join(str_path_proj, str_nm + '.seq')
    str_path_fesrep = os.path.join(opts.tdir, str_nm + '_fesrep.xml')
    str_proj = 'model'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_model = os.path.join(str_path_proj, str_nm + '.pkl')
    if opts.med == 'lasso':
        if opts.auto == False:
            print('Selecting sequence feature...')
            seqfs.LassoFeatureSelection(str_path_sg, str_path_lab, opts.ref, 
                                        str_path_seq, str_path_fesrep,
                                        str_path_model, opts.ups, opts.dws,
                                        opts.cv, opts.niter, opts.njob)
            print('Done')
        else:
            print('Detecting Optimal Feature Space...')
            seqfs.LassoRegionOptimizer(str_path_sg, str_path_lab, opts.ref,
                                       str_path_seq, str_path_fesrep,
                                       str_path_model, opts.ir,
                                       opts.ir+opts.rad, opts.step, opts.cv,
                                       opts.niter, opts.njob)
            print('Done')
    elif opts.med == 'logit':
        if opts.auto == False:
            print('Selecting sequence feature...')
            seqfs.LogitFeatureSelection(str_path_sg, str_path_lab, opts.ref, 
                                        str_path_seq, str_path_fesrep,
                                        str_path_model, opts.ups, opts.dws,
                                        opts.cv, opts.niter, opts.njob)
            print('Done')
        else:
            print('Detecting Optimal Feature Space...')
            seqfs.LogitRegionOptimizer(str_path_sg, str_path_lab, opts.ref,
                                       str_path_seq, str_path_fesrep,
                                       str_path_model, opts.ir,
                                       opts.ir+opts.rad, opts.step, opts.cv,
                                       opts.niter, opts.njob)
            print('Done')
        


