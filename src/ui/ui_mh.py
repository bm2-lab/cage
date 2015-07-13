import argparse

def ParseMh(p_mh):
    mh_inp = p_mh.add_argument_group('Input options')
    mh_inp.add_argument('-i', '--samind', dest='samind', required=True, help='sgRNA-Indel Table (required)')

    mh_otp = p_mh.add_argument_group('Output options')
    mh_otp.add_argument('-o', '--output-dir', dest='tdir', help='Output Directory, default = .', default='.')
    mh_par = p_mh.add_argument_group('Parameter options')
    mh_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    mh_par.add_argument('-t', '--cutoff', dest='cut', help='Cutoff of Reads (default = 0)', default=0, type=int)





