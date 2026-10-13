#call like this
#mamba install -c conda-forge compilers
#mamba install -c conda-forge cython numpy
#python setup.py build_ext --inplace
#python flow1_c_setup.py build_ext --inplace

#based on example here
#http://cython.readthedocs.io/en/latest/src/reference/compilation.html
import sys
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

#For Microsoft compiler, use '/openmp' instead of '-fopenmp'.
o_mp = '/openmp' if sys.platform == 'win32' else '-fopenmp'

extensions = [
    Extension("*", ["*.pyx"],
        include_dirs = [numpy.get_include()],
        extra_compile_args = [o_mp],
        extra_link_args=[o_mp],
        #define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
             ),
]

setup(
    name = "",
    ext_modules = cythonize(extensions),
)

