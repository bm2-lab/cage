import argparse
import os
import sys
from src import ui
        
def main():
    
    str_prog = 'CAGE'
    str_usage = 'python cage.py <command> [options]'
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

    p.add_argument('-v', '--version', action='version', version='%(prog)s 3.0.0(dev)')

    sp = p.add_subparsers(title='Command', metavar='')

    p_sg = sp.add_parser('sg',
                         description='Process sgRNA sequences into sgRNA information table',
                         usage='python cage.py sg [options]',
                         help='Processing sgRNA')
    ui.ParseSg(p_sg)
    
    p_prep = sp.add_parser('prep',
                           description='Process raw sequence data into sgRNA-Indel Table',
                           usage='python cage.py prep [options]',
                           help='Processing raw reads')

    ui.ParsePrep(p_prep)

    p_mh = sp.add_parser('mh',
                         description='Microhomology Detection',
                         usage='python cage.py mh [options]',
                         help='Microhomology Detection')
    ui.ParseMh(p_mh)
    
    p_ind = sp.add_parser('indel',
                          description='Indel Frameshifting Paradigm Analysis',
                          usage='python cage.py indel [options]',
                          help='Indel Frameshifting Paradigm Analysis')
    ui.ParseInd(p_ind)
    
    p_fs = sp.add_parser('fs',
                          description='Seq-Feature Selection and Model Prediction',
                          usage='python cage.py fs [options]',
                          help='Seq-Feature Selection and Model Prediction')
    ui.ParseFs(p_fs)

    p_mt = sp.add_parser('mt',
                         description='Seq-Feature Selection with Multi-task Group LASSO',
                         usage='python cage.py mt [options]',
                         help='Seq-Feature Selection with Multi-task Group LASSO')
    ui.ParseMt(p_mt)

    p_ev = sp.add_parser('eval',
                         description='sgRNA KO efficiency evaluation',
                         usage='python cage.py eval [options]',
                         help='sgRNA KO efficiency evaluation')
    ui.ParseEval(p_ev)
    
    p_vis = sp.add_parser('vis',
                          description='Visualization',
                          usage='python cage.py vis [options]',
                          help='Visualization')
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
        elif sys.argv[1] == 'fs':
            p_fs.print_help()
        elif sys.argv[1] == 'mt':
            p_mt.print_help()
        elif sys.argv[1] == 'eval':
            p_ev.print_help()
        elif sys.argv[1] == 'vis':
            p_vis.print_help()
        else:
            p.print_help()
        sys.exit(1)

    if 'BWADB' not in os.environ:
        p.error('$BWADB Not Exist. See README')
    
    if 'FASTADB' not in os.environ:
        p.error('$FASTADB Not Exist. See README')
    
    opts = p.parse_args()

    if not os.path.exists(opts.tdir):
            print('Target Directory not find, create one.')
            os.makedirs(opts.tdir)
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
    elif sys.argv[1] == 'fs':
        from src.core import fs
        fs.GeneralFs(opts)
    elif sys.argv[1] == 'mt':
        from src.core import mt
        mt.MtFs(opts)
    elif sys.argv[1] == 'eval':
        from src.core import eval
        eval.SgEval(opts)
    elif sys.argv[1] == 'vis':
        from src.core import vis
        vis.Visualize(opts)
        

if __name__ == '__main__':
    main()
