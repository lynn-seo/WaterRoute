from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("flow1_c",
                         ["flow1_c.pyx"],
                         include_dirs = [numpy.get_include()]
                        )]

setup(
  name = 'flow1_c',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
