import os
from src.core import prep
from src.core import sg
from src.core import seqfs
from sgscanner import *

def SgEval(opts):
    str_nm = os.path.basename(opts.tdir)
    str_path_mdl = opts.sfunc
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    str_path_sg = None
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
        
    if  opts.sg is None:
        dict_drct = {'two-sided': 'b', 'pos': 'f', 'neg': 'r'}
        str_path_sgfq = os.path.join(str_path_proj, str_nm + '.sgfq')
        gn_sg = ExtractSg(opts.ref, opts.chrom, opts.beg, opts.end, dict_drct[opts.drct])
        FaToFq(gn_sg, str_path_sgfq)
        
        str_path_sgpsam = os.path.join(str_path_proj, str_nm + '.sgpsam')
        str_path_sgsam = os.path.join(str_path_proj, str_nm + '.sgsam')
        str_path_sg = os.path.join(opts.tdir, str_nm + '.sg')
        print('Mapping sgRNA seq to ref genome with Bwa...')
        prep.CallBWA(str_path_sgfq, '', opts.ref, str_path_sgpsam, False, opts.thrd)
        prep.FilterSam(str_path_sgpsam, str_path_sgsam, False)
        print('Done')
        print('Processing sgsam...')
        sg.OrganizeSgsam(str_path_sgsam, str_path_sg)
        print('Done')
    else:
        str_path_sg = opts.sg
        
    str_path_seq = os.path.join(str_path_proj, str_nm + '.seq')
    str_path_st = os.path.join(opts.tdir, str_nm + '.st')
    print('sgRNA Evaluating...')
    seqfs.FeatureEval(str_path_sg, str_path_mdl, opts.ref,
                      str_path_seq, str_path_st)
    print('Done')
        


