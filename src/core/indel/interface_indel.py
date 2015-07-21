import os
from indel_integrator import *

def AnalyzeIndel(opts):
    str_nm = os.path.basename(os.path.splitext(opts.samind)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_samind = opts.samind
    str_path_iost = os.path.join(opts.tdir, str_nm + '.iost')
    str_path_iosg = os.path.join(str_path_proj, str_nm + '.iosg')
    print('Integrating Indel information...')
    AnalyzeSamind(str_path_samind, opts.sg, str_path_iost, str_path_iosg, opts.cut)
    print('Done')


