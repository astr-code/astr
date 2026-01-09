from numpy.distutils.core import setup, Extension

ext = Extension(
    name='p4pastr._pastr',       # compiled module inside package
    sources=['src/pastr_constdef.F90',
             'p4pastr/kernels.f90'], # Fortran source
             )

setup(
    name='p4pastr',
    version='0.1.0',
    description='Python analysis and driver package',
    author='Jian Fang',
    packages=['p4pastr'],
    ext_modules=[ext],
    python_requires='>=3.9',
)
