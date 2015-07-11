import os
from src.core import prep
from sgprocessor import *

def ProcessSg(p, opts):
    if opts.anno == True:
        if 'BEDDB' not in os.environ:
            p.error('$BEDDB Not Exist. See README')
    str_path_sgfq = opts.sg
    str_nm = os.path.basename(os.path.splitext(opts.sg)[0])
    str_proj = 'aux'
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.makedirs(str_path_proj)
    str_path_sgpsam = os.path.join(str_path_proj, str_nm + '.sgpsam')
    str_path_sgsam = os.path.join(str_path_proj, str_nm + '.sgsam')
    str_path_sg = os.path.join(opts.tdir, str_nm + '.sg')
    print('Mapping sgRNA seq to ref genome with Bwa...')
    prep.CallBWA(str_path_sgfq, '', opts.ref, str_path_sgpsam, False, opts.thrd)
    prep.FilterSam(str_path_sgpsam, str_path_sgsam, False)
    print('Done')
    print('Processing sgsam...')
    OrganizeSgsam(str_path_sgsam, str_path_sg)
    print('Done')
    if opts.anno == True:
        str_path_sgbed = os.path.join(str_path_proj, str_nm + '.sgbed')
        str_path_sgmap = os.path.join(str_path_proj, str_nm + '.sgmap')
        str_path_sga = os.path.join(opts.tdir, str_nm + '.sga')
        print('Annotating sgRNA...')
        int_status = AnnotateSg(str_path_sgsam, opts.ref, str_path_sgbed, str_path_sgmap)
        if int_status == 1:
            print('Annotated with RefSeq')
        elif int_status ==2:
            print('Annotated with RefSeq and UCSC Gene')
        elif int_status ==3:
            print('Annotated with RefSeq, UCSC Gene and GENCODE')
        elif int_status == 4:
            print('Annotated with RefSeq and UCSC Gene')
            print('Warning: Some are marked with None')
        elif int_status == 5:
            print('Annotated with RefSeq, UCSC Gene and GENCODE')
            print('Warning: Some are marked with None')
        print('Done')
        print('Merging sg and sgmap...')
        MergeSg(str_path_sg, str_path_sgmap, str_path_sga)
        print('Done')
