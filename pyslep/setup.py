from distutils.core import setup, Extension
from Cython.Build import cythonize
from numpy import get_include

ext_main = Extension('mtLeastR', ['mtLeastR/mtLeastR.pyx', 'mtLeastR/ep21R.c', 'mtLeastR/eppMatrix.c', 'mtLeastR/epph.c'], include_dirs=['.', 'mtLeastR', get_include()])

setup(ext_modules=cythonize([ext_main]))
