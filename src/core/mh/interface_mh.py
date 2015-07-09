import os
from mhprocessor import *

def DetectMh(opts):
    str_nm = os.path.basename(os.path.splitext(opts.samind)[0])
    str_path_proj = opts.tdir
    str_path_samind = opts.samind
    str_path_mh = os.path.join(opts.tdir, str_nm + '.mh')
    str_path_mnst = os.path.join(opts.tdir, str_nm + '.mnst')
    print('Detecting microhomology...')
    ProcessMh(str_path_samind, opts.ref, str_path_mh)
    print('Done')
    print('Analyzing microhomology...')
    AnalyzeMh(str_path_mh, str_path_mnst, opts.cut)
    print('Done')


