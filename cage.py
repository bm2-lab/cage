import argparse
import os
import sys
import subprocess

from src import prep
from src import mh
from src import indel
from src import vis

def Preprocess(opts):
    bln_pair = True if opts.rev else False
    opts.rev = opts.rev if bln_pair else ''
    bln_long = True if opts.length == 'long' else False
    str_proj = os.path.basename(os.path.splitext(opts.fwd)[0])[:-2] if bln_pair else os.path.basename(os.path.splitext(opts.fwd)[0])
    str_path_proj = os.path.join(opts.tdir, str_proj)
    if not os.path.exists(str_path_proj):
        os.mkdir(str_path_proj)
    str_path_sam = os.path.join(str_path_proj, 'sam')
    if not os.path.exists(str_path_sam):
        os.mkdir(str_path_sam)
    str_path_psam = os.path.join(str_path_sam, str_proj + '.psam')
    str_path_sam = os.path.join(str_path_sam, str_proj + '.sam')

    str_fsg = os.path.basename(os.path.splitext(opts.sg)[0])
    str_path_sg = os.path.join(str_path_proj, 'sg')
    if not os.path.exists(str_path_sg):
        os.mkdir(str_path_sg)
    str_path_sgfq = opts.sg
    str_path_sgpsam = os.path.join(str_path_sg, str_fsg + '.sgpsam')
    str_path_sgsam = os.path.join(str_path_sg, str_fsg + '.sgsam')
    str_path_sgbed = os.path.join(str_path_sg, str_fsg + '.sgbed')
    str_path_sgmap = os.path.join(str_path_sg, str_fsg + '.sgmap')
    str_path_psg = os.path.join(str_path_sg, str_fsg + '.psg')
    str_path_sg = os.path.join(opts.tdir, str_fsg + '.sg')

    str_path_samind = os.path.join(opts.tdir, str_proj + '.samind')

    print('Mapping sgRNA seq to ref genome with Bwa...')
    prep.CallBWA(str_path_sgfq, '', opts.ref, str_path_sgpsam, False, opts.thrd)
    prep.FilterSam(str_path_sgpsam, str_path_sgsam, False)
    print('Done')

    print('Annotating sgRNA...')
    int_status = prep.AnnotateSg(str_path_sgsam, opts.ref, str_path_sgbed, str_path_sgmap)
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

    print('Processing sgsam...')
    prep.OrganizeSgsam(str_path_sgsam, str_path_psg)
    print('Done')

    print('Merging psg and sgmap...')
    prep.MergeSg(str_path_psg, str_path_sgmap, str_path_sg)
    print('Done')

    print('Mapping target reads with Bwa...')
    prep.CallBWA(opts.fwd, opts.rev, opts.ref, str_path_psam, bln_long, opts.thrd)
    prep.FilterSam(str_path_psam, str_path_sam, bln_pair)
    print('Done')

    print('Integrating sam with sg...')
    prep.ProcessSam(str_path_sam, str_path_sg, str_path_samind, bln_pair,
               int_idloffset=20, int_idroffset=20, int_loffset=0, int_roffset=0)
    print('Done')

    

def DetectMh(opts):
    str_nm = os.path.basename(os.path.splitext(opts.samind)[0])
    str_path_proj = opts.tdir
    str_path_samind = opts.samind
    str_path_mh = os.path.join(opts.tdir, str_nm + '.mh')
    str_path_mnst = os.path.join(opts.tdir, str_nm + '.mnst')
    print('Detecting microhomology...')
    mh.ProcessMh(str_path_samind, opts.ref, str_path_mh)
    print('Done')
    print('Analyzing microhomology...')
    mh.AnalyzeMh(str_path_mh, str_path_mnst)
    print('Done')

def AnalyzeIndel(opts):
    str_nm = os.path.basename(os.path.splitext(opts.samind)[0])
    str_path_proj = opts.tdir
    str_path_samind = opts.samind
    str_path_iost = os.path.join(opts.tdir, str_nm + '.iost')
    str_path_seq = os.path.join(opts.tdir, str_nm + '.seq')
    str_path_fesrep = os.path.join(opts.tdir, str_nm + '_fesrep.xml')
    print('Integrating Indel information...')
    indel.AnalyzeSamind(str_path_samind, str_path_iost, opts.cut)
    print('Done')
    print('Extracting sequence feature...')
    indel.ExtractSeqFeature(opts.sg, opts.ref, str_path_seq, opts.ups, opts.dws)
    print('Done')
    print('Selecting sequence feature...')
    indel.LassoFeatureSelection(str_path_seq, str_path_iost, str_path_fesrep, opts.cv, opts.niter)
    print('Done')

def Visualize(opts):
    str_path_proj = opts.tdir
    if opts.fes != None:
        str_nm = os.path.basename(os.path.splitext(opts.fes)[0])
        str_path_fes = opts.fes
        str_path_tex = os.path.join(opts.tdir, str_nm + '.tex')
        print('Visualizing selected feature...')
        vis.VisualizeFeature(str_path_fes, str_path_tex, opts.amp)
        str_latex = r'pdflatex -output-directory=%s %s'% (opts.tdir, str_path_tex)
        subprocess.call(str_latex, shell=True, executable='/bin/bash')
        print('Done')

def main():
    
    str_prog = 'CAGE'
    str_usage = 'cage <command> [options]'
    str_desc = r'''
CRISPR KO Analysis based on Genomic Editing data
------------------------------------------------
A CRISPR-cas9 based Genome Editing data analysis pipeline, 
for the analysis of indels and microhomology patterns from 
CRISPR-Cas9 Knock-Out NGS data.
'''
    str_epilog = r'Note: Environment Variable $BWADB $BEDDB $FASTADB Must Exist.'


    p = argparse.ArgumentParser(prog=str_prog,
                                 usage=str_usage,
                                 description=str_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=str_epilog)

    p.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.2(dev)')

    sp = p.add_subparsers(title='Command', metavar='')

    p_prep = sp.add_parser('prep',
                           description='Process raw sequence data into sgRNA-Indel Table',
                           usage='cage prep [options]',
                           help='Preprocessing')

    

    prep_inp = p_prep.add_argument_group('Input options')
    prep_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Fastq File (required)')
    prep_inp.add_argument('-f', '--forward', dest='fwd', required=True, help='Forward Fastq File (required)')
    prep_inp.add_argument('-r', '--reverse', dest='rev', help='Reverse Fastq File')

    prep_otp = p_prep.add_argument_group('Output options')
    prep_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')

    prep_par = p_prep.add_argument_group('Parameter options')
    prep_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    prep_par.add_argument('-l', '--length', dest='length', help='Read Length (long > 70bp, short <= 70bp. default = long)', choices=('long', 'short'), default='long')
    prep_par.add_argument('-t', '--thread', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int)

    p_mh = sp.add_parser('mh',
                         description='Microhomology Detection',
                         usage='cage mh [options]',
                         help='Microhomology Detection')
    mh_inp = p_mh.add_argument_group('Input options')
    mh_inp.add_argument('-i', '--samind', dest='samind', required=True, help='sgRNA-Indel Table (required)')

    mh_otp = p_mh.add_argument_group('Output options')
    mh_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')
    mh_par = p_mh.add_argument_group('Parameter options')
    mh_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    
    p_ind = sp.add_parser('indel',
                          description='Lasso for Indel Frameshifting Paradigm',
                          usage='cage indel [options]',
                          help='Indel Feature Selection')
    ind_inp = p_ind.add_argument_group('Input options')
    ind_inp.add_argument('-i', '--samind', dest='samind', required=True, help='sgRNA-Indel Table (required)')
    ind_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')

    ind_otp = p_ind.add_argument_group('Output options')
    ind_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')
    
    ind_par = p_ind.add_argument_group('Parameter options')
    ind_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    ind_par.add_argument('-t', '--cutoff', dest='cut', help='Cutoff of Reads (default = 0)', default=0, type=int)
    ind_par.add_argument('-u', '--ups', dest='ups', help='Upstream Region Length (default = 35)', default=35, type=int)
    ind_par.add_argument('-w', '--dws', dest='dws', help='Downstream Region Length (default = 32)', default=32, type=int)

    ind_las = p_ind.add_argument_group('Feature Selection options')
    ind_las.add_argument('-c', '--cv', dest='cv', help='Folds for Cross Validation (default = 10)', default=10, type=int)
    ind_las.add_argument('-n', '--niter', dest='niter', help='Iteration Times for Cross Validation (default = 3000)', default=3000, type=int)

    p_vis = sp.add_parser('vis',
                          description='Result Visualization',
                          usage='cage vis [options]',
                          help='Result Visualization')
    vis_fes = p_vis.add_argument_group('Feature Visualization options')
    vis_fes.add_argument('-f', '--fes', dest='fes', help='Feature Report File')
    vis_fes.add_argument('-a', '--amp', dest='amp', help='Axis Amplifier (default = 1.00)', default=1.00, type=float)
    vis_otp = p_vis.add_argument_group('Output options')
    vis_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')
    
    

    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'prep':
            p_prep.print_help()
        elif sys.argv[1] == 'mh':
            p_mh.print_help()
        elif sys.argv[1] == 'indel':
            p_ind.print_help()
        elif sys.argv[1] == 'vis':
            p_vis.print_help()
        sys.exit(1)

    if 'BWADB' not in os.environ:
        p.error('$BWADB Not Exist. See --help')
    if 'BEDDB' not in os.environ:
        p.error('$BEDDB Not Exist. See --help')
    if 'FASTADB' not in os.environ:
        p.error('$FASTADB Not Exist. See --help')
    
    opts = p.parse_args()

    if not os.path.exists(opts.tdir):
            print('Target Directory not find, create one.')
            os.mkdir(opts.tdir)
    if sys.argv[1] == 'prep':
        Preprocess(opts)
    elif sys.argv[1] == 'mh':
        DetectMh(opts)
    elif sys.argv[1] == 'indel':
        AnalyzeIndel(opts)
    elif sys.argv[1] == 'vis':
        Visualize(opts)
        

if __name__ == '__main__':
    main()
