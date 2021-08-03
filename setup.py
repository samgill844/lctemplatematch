from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy
ext_modules = [
    Extension(
        "cybinary",
        ["source/kepler.pyx"],#, "source/rv_code.pyx"],
    )
]

setup(
    name='cybinary',
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-O3','-fopenmp'],
    extra_link_args=['-fopenmp'],
    scripts=['Utils/tls']
)