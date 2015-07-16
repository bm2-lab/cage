import argparse

def ParseEval(p_ev):
    ev_inp = p_ev.add_argument_group('Input options')
    ev_inp.add_argument('-c', '--chrom', dest='chrom', help='Target Chromosome')
    ev_inp.add_argument('-b', '--beg', dest='beg', help='Start Coordinate', type=int)
    ev_inp.add_argument('-e', '--end', dest='end', help='End Coordinate', type=int)
    ev_inp.add_argument('-s', '--sg', dest='sg', help='sgRNA Information File')
    ev_inp.add_argument('-f', '--score-function', dest='sfunc', required=True, help='Score Function (required)')

    ev_otp = p_ev.add_argument_group('Output options')
    ev_otp.add_argument('-o', '--output-dir', dest='tdir', help='Output Directory, default = .', default='.')

    ev_par = p_ev.add_argument_group('Parameter options')
    ev_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    ev_par.add_argument('-d', '--direction', dest='drct', help='Scan Direction (default = two-sided)', default='two-sided', choices=('two-sided', 'pos', 'neg'))
    ev_par.add_argument('-t', '--thread', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int)
    



