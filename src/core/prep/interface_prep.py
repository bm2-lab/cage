import os
from fqmapper import *
from samprocessor import *

def Preprocess(opts):
    bln_pair = True if opts.rev else False
    opts.rev = opts.rev if bln_pair else ''
    bln_long = True if opts.length == 'long' else False
    str_proj = os.path.basename(os.path.splitext(opts.fwd)[0])[:-2] if bln_pair else os.path.basename(os.path.splitext(opts.fwd)[0])
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
        
    str_path_sg = opts.sg
    str_path_sam = os.path.join(str_path_proj, 'sam')
    if not os.path.exists(str_path_sam):
        os.makedirs(str_path_sam)
    str_path_psam = os.path.join(str_path_sam, str_proj + '.psam')
    str_path_sam = os.path.join(str_path_sam, str_proj + '.sam')
    str_path_samind = os.path.join(opts.tdir, str_proj + '.samind')

    print('Mapping raw reads with Bwa...')
    CallBWA(opts.fwd, opts.rev, opts.ref, str_path_psam, bln_long, opts.thrd)
    FilterSam(str_path_psam, str_path_sam, bln_pair)
    print('Done')

    print('Integrating sam with sg...')
    ProcessSam(str_path_sam, str_path_sg, str_path_samind, bln_pair,
               int_idloffset=20, int_idroffset=20, int_loffset=0, int_roffset=0)
    print('Done')


