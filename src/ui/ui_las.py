import argparse

def ParseLas(p_las):
    las_inp = p_las.add_argument_group('Input options')
    las_inp.add_argument('-i', '--input-label', dest='label', required=True, help='Input Label Dataset (required)')
    las_inp.add_argument('-s', '--sg', dest='sg', required=True, help='sgRNA Information File (required)')

    las_otp = p_las.add_argument_group('Output options')
    las_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')
    
    las_par = p_las.add_argument_group('Parameter options')
    las_par.add_argument('-g', '--reference', dest='ref', required=True, help='Reference Genome (required)')
    las_par.add_argument('-a', '--auto', dest='auto', action='store_true', help='Auto Detection for Sequence Region')
    las_par.add_argument('--init-radius', dest='ir', help='Detection Region Init Radius (default = 0)', default=0, type=int)
    las_par.add_argument('-r', '--radius', dest='rad', help='Detection Region Radius (default = 200)', default=200, type=int)
    las_par.add_argument('--step', dest='step', help='Detection Step (default = 5)', default=5, type=int)
    las_par.add_argument('-u', '--ups', dest='ups', help='Upstream Region Length (default = 30)', default=30, type=int)
    las_par.add_argument('-w', '--dws', dest='dws', help='Downstream Region Length (default = 27)', default=27, type=int)

    las_las = p_las.add_argument_group('Feature Selection options')
    las_las.add_argument('-c', '--cv', dest='cv', help='Folds for Cross Validation (default = 5)', default=5, type=int)
    las_las.add_argument('-n', '--niter', dest='niter', help='Iteration Times for Cross Validation (default = 1000)', default=1000, type=int)





