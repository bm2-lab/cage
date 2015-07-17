import argparse

def ParseSg(p_sg):
    sg_inp = p_sg.add_argument_group('Input options')
    sg_inp.add_argument('-s', dest='sg', required=True, help='sgRNA Fastq File (required)', metavar='<sgRNA.fq>')
    sg_otp = p_sg.add_argument_group('Output options')
    sg_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')

    sg_par = p_sg.add_argument_group('Parameter options')
    sg_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    sg_par.add_argument('-a', dest = 'anno', action='store_true', help='Annotate sgRNA sequences (optional, if annotation is needed)')
    sg_par.add_argument('-t', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int, metavar='<int>')




