import argparse

def ParsePrep(p_prep):
    prep_inp = p_prep.add_argument_group('Input options')
    prep_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Fastq File (required)')
    prep_inp.add_argument('-f', '--forward', dest='fwd', required=True, help='Forward Fastq File (required)')
    prep_inp.add_argument('-r', '--reverse', dest='rev', help='Reverse Fastq File')

    prep_otp = p_prep.add_argument_group('Output options')
    prep_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')

    prep_par = p_prep.add_argument_group('Parameter options')
    prep_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    prep_par.add_argument('-l', '--length', dest='length', help='Read Length (long > 70bp, short <= 70bp. default = long)', choices=('long', 'short'), default='long')
    prep_par.add_argument('-t', '--thread', dest='thrd', help='Threads for BWA, default = 1', default=1, type=int)




