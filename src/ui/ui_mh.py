import argparse

def ParseMh(p_mh):
    mh_inp = p_mh.add_argument_group('Input options')
    mh_inp.add_argument('-i', dest='samind', required=True, help='sgRNA-Indel Table (required)', metavar='<samind file>')

    mh_otp = p_mh.add_argument_group('Output options')
    mh_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')
    mh_par = p_mh.add_argument_group('Parameter options')
    mh_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    mh_par.add_argument('-t', dest='cut', help='Read-Count Cutoff (default = 0)', default=0, type=int, metavar='<int>')





