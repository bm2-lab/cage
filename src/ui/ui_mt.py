import argparse

def ParseMt(p_mt):
    mt_inp = p_mt.add_argument_group('Input options')
    mt_inp.add_argument('-i', nargs='*', dest='label', required=True, help='Label Dataset List (required)', metavar='<label file>')
    mt_inp.add_argument('-s', nargs='*', dest='sg', required=True, help='sgRNA Information Table List (required)', metavar='<sg file>')

    mt_otp = p_mt.add_argument_group('Output options')
    mt_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')
    
    mt_par = p_mt.add_argument_group('Parameter options')
    mt_par.add_argument('-g', nargs='*', dest='ref', required=True, help='Reference Genome List (required)', metavar='<ref genome>')
    mt_par.add_argument('-u', dest='ups', help='Upstream Region Length (default = 30)', default=30, type=int, metavar='<int>')
    mt_par.add_argument('-w', dest='dws', help='Downstream Region Length (default = 27)', default=27, type=int, metavar='<int>')
    mt_par.add_argument('-d', dest='lmd', help='Selection Strictness (between 0 and 1, default = 0.5)', default=0.5, type=float, metavar='<float>')






