import argparse

def ParseInd(p_ind):
    ind_inp = p_ind.add_argument_group('Input options')
    ind_inp.add_argument('-i', dest='samind', required=True, help='sgRNA-Indel Table (required)', metavar='<samind file>')
    ind_inp.add_argument('-s', dest='sg', required=True, help='sgRNA Information Table (required)', metavar='<sg file>')

    ind_otp = p_ind.add_argument_group('Output options')
    ind_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')
    
    ind_par = p_ind.add_argument_group('Parameter options')
    ind_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    ind_par.add_argument('-t', dest='cut', help='Read-Count Cutoff (default = 0)', default=0, type=int, metavar='<int>')





