import argparse

def ParseFs(p_fs):
    fs_inp = p_fs.add_argument_group('Input options')
    fs_inp.add_argument('-i', dest='label', required=True, help='Label Dataset (required)', metavar='<label file>')
    fs_inp.add_argument('-s', dest='sg', required=True, help='sgRNA Information Table (required)', metavar='<sg file>')

    fs_otp = p_fs.add_argument_group('Output options')
    fs_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')
    
    fs_par = p_fs.add_argument_group('Parameter options')
    fs_par.add_argument('-g', dest='ref', required=True, help='Reference Genome (required)', metavar='<ref genome>')
    fs_par.add_argument('-a', dest='auto', action='store_true', help='Auto Detection for Optimal Sequence Region')
    fs_par.add_argument('--init-radius', dest='ir', help='Detection Region Init Radius (default = 0)', default=0, type=int, metavar='<int>')
    fs_par.add_argument('-r', dest='rad', help='Detection Region Radius (default = 200)', default=200, type=int, metavar='<int>')
    fs_par.add_argument('--step', dest='step', help='Detection Step (default = 5)', default=5, type=int, metavar='<int>')
    fs_par.add_argument('-u', dest='ups', help='Upstream Region Length (default = 30)', default=30, type=int, metavar='<int>')
    fs_par.add_argument('-w', dest='dws', help='Downstream Region Length (default = 27)', default=27, type=int, metavar='<int>')

    fs_fs = p_fs.add_argument_group('Feature Selection options')
    fs_fs.add_argument('-m', dest='med', required=True, help='Method Selection (required)', choices=('lasso', 'logit'), metavar='<lasso|logit>')
    fs_fs.add_argument('-c', dest='cv', help='Folds for Cross Validation (default = 5)', default=5, type=int, metavar='<int>')
    fs_fs.add_argument('-n', dest='niter', help='Iteration Times for Cross Validation (default = 1000)', default=1000, type=int, metavar='<int>')
    fs_fs.add_argument('-j', dest='njob', help='Number of CPU cores used (default = 1; -1 means all cores)', default=1, type=int, metavar='<int>')





