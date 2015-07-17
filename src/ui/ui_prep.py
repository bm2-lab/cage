import argparse

def ParsePrep(p_prep):
    prep_inp = p_prep.add_argument_group('Input options')
    prep_inp.add_argument('-s', dest='sg', required=True, help='sgRNA Information Table (required)', metavar='<sg file>')
    prep_inp.add_argument('-f', dest='fwd', required=True, help='Forward Fastq File (required)', metavar='<reads_1.fq>')
    prep_inp.add_argument('-r', dest='rev', help='Reverse Fastq File', metavar='<reads_2.fq>')

    prep_otp = p_prep.add_argument_group('Output options')
    prep_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')

    prep_par = p_prep.add_argument_group('Parameter options')
    prep_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    prep_par.add_argument('-l', '--length', dest='length', help='Read Length (long > 70bp, short <= 70bp. default = long)', choices=('long', 'short'), default='long', metavar='<long|short>')
    prep_par.add_argument('-t', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int, metavar='<int>')




