import argparse

def ParseFs(p_fs):
    fs_inp = p_fs.add_argument_group('Input options')
    fs_inp.add_argument('-i', '--input-label', dest='label', required=True, help='Input Label Dataset (required)')
    fs_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')

    fs_otp = p_fs.add_argument_group('Output options')
    fs_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')
    
    fs_par = p_fs.add_argument_group('Parameter options')
    fs_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    fs_par.add_argument('-a', '--auto', dest='auto', action='store_true', help='Auto Detection for Optimal Sequence Region')
    fs_par.add_argument('--init-radius', dest='ir', help='Detection Region Init Radius (default = 0)', default=0, type=int)
    fs_par.add_argument('-r', '--radius', dest='rad', help='Detection Region Radius (default = 200)', default=200, type=int)
    fs_par.add_argument('--step', dest='step', help='Detection Step (default = 5)', default=5, type=int)
    fs_par.add_argument('-u', '--ups', dest='ups', help='Upstream Region Length (default = 30)', default=30, type=int)
    fs_par.add_argument('-w', '--dws', dest='dws', help='Downstream Region Length (default = 27)', default=27, type=int)

    fs_fs = p_fs.add_argument_group('Feature Selection options')
    fs_fs.add_argument('-m', '--method', dest='med', required=True, help='Method Selection (required)', choices=('lasso', 'logit'))
    fs_fs.add_argument('-c', '--cv', dest='cv', help='Folds for Cross Validation (default = 5)', default=5, type=int)
    fs_fs.add_argument('-n', '--niter', dest='niter', help='Iteration Times for Cross Validation (default = 1000)', default=1000, type=int)
    fs_fs.add_argument('-j', '--jobs', dest='njob', help='Number of CPU cores used (default = 1; -1 means all cores)', default=1, type=int)





