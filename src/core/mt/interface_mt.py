import os
from src.core import seqfs

def MtFs(opts):
    func_bnm = lambda x: os.path.basename(os.path.splitext(x)[0])
    func_nm = lambda x, y: '%s_%s'% (func_bnm(x), func_bnm(y))
    str_nm = reduce(func_nm, opts.label)
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    lst_path_lab = opts.label
    lst_path_sg = opts.sg
    lst_path_ref = opts.ref
    lst_seq = ['%s.seq'%i for i in map(func_bnm, opts.sg)]
    lst_path_seq = [os.path.join(str_path_proj, i) for i in lst_seq]
    str_path_fesrep = os.path.join(opts.tdir, str_nm + '_fesrep.xml')
    print('Selecting sequence feature...')
    seqfs.MtLassoFeatureSelection(lst_path_sg, lst_path_lab, lst_path_ref,
                                  lst_path_seq, str_path_fesrep,
                                  opts.ups, opts.dws, opts.lmd)
    print('Done')
