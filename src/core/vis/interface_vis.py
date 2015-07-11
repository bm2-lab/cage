import os
import subprocess
from feature_visualizer import *

def Visualize(opts):
    str_path_proj = opts.tdir
    if opts.fes != None:
        str_nm = os.path.basename(os.path.splitext(opts.fes)[0])
        str_path_fes = opts.fes
        str_path_tex = os.path.join(opts.tdir, str_nm + '.tex')
        print('Visualizing selected feature...')
        VisualizeFeature(str_path_fes, str_path_tex, opts.amp)
        str_latex = r'pdflatex -output-directory=%s %s'% (opts.tdir, str_path_tex)
        subprocess.call(str_latex, shell=True, executable='/bin/bash')
        print('Done')



