import os
from src.core import seqfs

def SgEval(opts):
    str_path_sg = opts.sg
    str_path_mdl = opts.sfunc
    str_nm = os.path.basename(os.path.splitext(opts.sg)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_seq = os.path.join(str_path_proj, str_nm + '.seq')
    str_path_st = os.path.join(opts.tdir, str_nm + '.st')
    print('sgRNA Evaluating...')
    seqfs.FeatureEval(str_path_sg, str_path_mdl, opts.ref,
                      str_path_seq, str_path_st)
    print('Done')
        


