import argparse

def ParseVis(p_vis):
    vis_fes = p_vis.add_argument_group('Feature Visualization options')
    vis_fes.add_argument('-f', '--fes', dest='fes', help='Feature Report File')
    vis_fes.add_argument('-a', '--amp', dest='amp', help='Axis Amplifier (default = 1.00)', default=1.00, type=float)
    vis_otp = p_vis.add_argument_group('Output options')
    vis_otp.add_argument('-d', '--target-dir', dest='tdir', help='Target Directory, default = .', default='.')





