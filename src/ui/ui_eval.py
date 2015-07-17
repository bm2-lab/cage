import argparse

def ParseEval(p_ev):
    ev_inp = p_ev.add_argument_group('Input options')
    ev_inp.add_argument('-c', dest='chrom', help='Target Chromosome', metavar='<chromosome>')
    ev_inp.add_argument('-b', dest='beg', help='Start Coordinate', type=int, metavar='<int>')
    ev_inp.add_argument('-e', dest='end', help='End Coordinate', type=int, metavar='<int>')
    ev_inp.add_argument('-s', dest='sg', help='sgRNA Information Table', metavar='<sg file>')
    ev_inp.add_argument('-f', dest='sfunc', required=True, help='Score Function (required)', metavar='<score function file>')

    ev_otp = p_ev.add_argument_group('Output options')
    ev_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')

    ev_par = p_ev.add_argument_group('Parameter options')
    ev_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    ev_par.add_argument('-d', dest='drct', help='Scan Direction (default = two-sided)', default='two-sided', choices=('two-sided', 'pos', 'neg'), metavar='<two-sided|pos|neg>')
    ev_par.add_argument('-t', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int, metavar='<int>')
    



