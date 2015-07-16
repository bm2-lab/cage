import argparse

def ParseInd(p_ind):
    ind_inp = p_ind.add_argument_group('Input options')
    ind_inp.add_argument('-i', '--samind', dest='samind', required=True, help='sgRNA-Indel Table (required)')
    ind_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')

    ind_otp = p_ind.add_argument_group('Output options')
    ind_otp.add_argument('-o', '--output-dir', dest='tdir', help='Output Directory, default = .', default='.')
    
    ind_par = p_ind.add_argument_group('Parameter options')
    ind_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    ind_par.add_argument('-t', '--cutoff', dest='cut', help='Read-Count Cutoff (default = 0)', default=0, type=int)





