import os
import glob
from setuptools import setup, find_packages, Extension

with open("README.md", encoding="utf-8") as fname:
    README = fname.read()


version = 'x.y.z'
if os.path.exists('VERSION'):
    version = open('VERSION').read().strip()

USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("homopolymer_compression", ["homopolymer_compression"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='tiptoft',
    version=version,
    description='tiptoft: predict which plasmid should be present from uncorrected long read data',
    long_description=README,
    packages = find_packages(),
    author='Andrew J. Page',
    author_email='andrew.page@quadram.ac.uk',
    url='https://github.com/andrewjpage/tiptoft',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'biopython >= 1.68',
        'pyfastaq >= 3.12.0',
        'cython'
    ],
    ext_modules = extensions,
    package_data={'tiptoft': ['data/*']},
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
