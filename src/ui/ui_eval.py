import argparse

def ParseEval(p_ev):
    ev_inp = p_ev.add_argument_group('Input options')
    ev_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')
    ev_inp.add_argument('-f', '--score-function', dest='sfunc', required=True, help='Score Function (required)')

    ev_otp = p_ev.add_argument_group('Output options')
    ev_otp.add_argument('-o', '--output-dir', dest='tdir', help='Output Directory, default = .', default='.')

    ev_par = p_ev.add_argument_group('Parameter options')
    ev_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')



