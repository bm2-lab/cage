import argparse

def ParseSg(p_sg):
    sg_inp = p_sg.add_argument_group('Input options')
    sg_inp.add_argument('-s', '--sg', required=True, help='sgRNA Fastq File (required)')
    sg_otp = p_sg.add_argument_group('Output options')
    sg_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')

    sg_par = p_sg.add_argument_group('Parameter options')
    sg_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    sg_par.add_argument('-a', '--annotate', dest = 'anno', action='store_true', help='Annotate sgRNA sequences (optional, if annotation is needed)')
    sg_par.add_argument('-t', '--thread', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int)




