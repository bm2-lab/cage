import argparse
import os
import sys
from src import ui
        
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
    p = argparse.ArgumentParser(prog=str_prog,
                                 usage=str_usage,
                                 description=str_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.2(dev)')

    sp = p.add_subparsers(title='Command', metavar='')

    p_sg = sp.add_parser('sg',
                         description='Process sgRNA sequences into sgRNA information table',
                         usage='cage sg [options]',
                         help='Processing sgRNA')
    ui.ParseSg(p_sg)
    
    p_prep = sp.add_parser('prep',
                           description='Process raw sequence data into sgRNA-Indel Table',
                           usage='cage prep [options]',
                           help='Processing raw reads')

    ui.ParsePrep(p_prep)

    p_mh = sp.add_parser('mh',
                         description='Microhomology Detection',
                         usage='cage mh [options]',
                         help='Microhomology Detection')
    ui.ParseMh(p_mh)
    
    p_ind = sp.add_parser('indel',
                          description='Lasso for Indel Frameshifting Paradigm',
                          usage='cage indel [options]',
                          help='Indel Feature Selection')
    ui.ParseInd(p_ind)

    p_vis = sp.add_parser('vis',
                          description='Result Visualization',
                          usage='cage vis [options]',
                          help='Result Visualization')
    ui.ParseVis(p_vis)

    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'sg':
            p_sg.print_help()
        elif sys.argv[1] == 'prep':
            p_prep.print_help()
        elif sys.argv[1] == 'mh':
            p_mh.print_help()
        elif sys.argv[1] == 'indel':
            p_ind.print_help()
        elif sys.argv[1] == 'vis':
            p_vis.print_help()
        sys.exit(1)

    if 'BWADB' not in os.environ:
        p.error('$BWADB Not Exist. See README')
    
    if 'FASTADB' not in os.environ:
        p.error('$FASTADB Not Exist. See README')
    
    opts = p.parse_args()

    if not os.path.exists(opts.tdir):
            print('Target Directory not find, create one.')
            os.mkdir(opts.tdir)
    if sys.argv[1] == 'sg':
        from src.core import sg
        sg.ProcessSg(p, opts)        
    elif sys.argv[1] == 'prep':
        from src.core import prep
        prep.Preprocess(opts)
    elif sys.argv[1] == 'mh':
        from src.core import mh
        mh.DetectMh(opts)
    elif sys.argv[1] == 'indel':
        from src.core import indel
        indel.AnalyzeIndel(opts)
    elif sys.argv[1] == 'vis':
        from src.core import vis
        vis.Visualize(opts)
        

if __name__ == '__main__':
    main()
