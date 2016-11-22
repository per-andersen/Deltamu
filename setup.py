from distutils.core import setup
from Cython.Build import cythonize

#Run with python setup.py build_ext --inplace
setup(
  name = 'scalar_inv_efuncs',
  ext_modules = cythonize("scalar_inv_efuncs.pyx"),
)