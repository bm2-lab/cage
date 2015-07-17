import argparse

def ParseVis(p_vis):
    vis_fes = p_vis.add_argument_group('Feature Visualization options')
    vis_fes.add_argument('-f', dest='fes', help='Feature Report File', metavar='<feature report>')
    vis_fes.add_argument('-a', dest='amp', help='Axis Amplifier (default = 1.00)', default=1.00, type=float, metavar='<float>')
    vis_otp = p_vis.add_argument_group('Output options')
    vis_otp.add_argument('-o', dest='tdir', help='Output Directory, default = .', default='.', metavar='<output directory>')





