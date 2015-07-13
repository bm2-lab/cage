import argparse

def ParseInd(p_ind):
    ind_inp = p_ind.add_argument_group('Input options')
    ind_inp.add_argument('-i', '--samind', dest='samind', required=True, help='sgRNA-Indel Table (required)')
    ind_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')

    ind_otp = p_ind.add_argument_group('Output options')
    ind_otp.add_argument('-o', '--output-dir', dest='tdir', help='Output Directory, default = .', default='.')
    
    ind_par = p_ind.add_argument_group('Parameter options')
    ind_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    ind_par.add_argument('-t', '--cutoff', dest='cut', help='Cutoff of Reads (default = 0)', default=0, type=int)
    ind_par.add_argument('-a', '--auto', dest='auto', action='store_true', help='Auto Detection for Sequence Region')
    ind_par.add_argument('--ir', dest='ir', help='Detection Region Init Radius (default = 0)', default=0, type=int)
    ind_par.add_argument('-r', '--radius', dest='rad', help='Detection Region Radius (default = 200)', default=200, type=int)
    ind_par.add_argument('--step', dest='step', help='Detection Step (default = 5)', default=5, type=int)
    ind_par.add_argument('-u', '--ups', dest='ups', help='Upstream Region Length (default = 30)', default=30, type=int)
    ind_par.add_argument('-w', '--dws', dest='dws', help='Downstream Region Length (default = 27)', default=27, type=int)

    ind_las = p_ind.add_argument_group('Feature Selection options')
    ind_las.add_argument('-c', '--cv', dest='cv', help='Folds for Cross Validation (default = 5)', default=5, type=int)
    ind_las.add_argument('-n', '--niter', dest='niter', help='Iteration Times for Cross Validation (default = 1000)', default=1000, type=int)
    ind_las.add_argument('-j', '--jobs', dest='njob', help='Number of CPU cores used (default = 1; -1 means all cores)', default=1, type=int)





