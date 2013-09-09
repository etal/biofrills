#!/usr/bin/env python

"""Bioinformatics utilities for molecular sequence analysis."""

from glob import glob
from os.path import join, isfile

setup_args = dict(
    name='biofrills',
    version='0.3.1',
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@gmail.com',
    url='http://github.com/etal/biofrills',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['biofrills', 'biofrills.stats'],
    #scripts=glob(join('scripts', '*'))
)

try:
    from setuptools import setup
    from setuptools.extension import Extension
    # Dependencies for easy_install and pip
    setup_args.update(
        install_requires=[
            'biopython >= 1.60',
            # 'Cython >= 0.15',
        ])
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

do_cython = False
try:
    from Cython.Distutils import build_ext
    if isfile('biofrills/cpairutils.pyx'):
        do_cython = True
except ImportError:
    pass

if do_cython:
    setup_args.update(
        cmdclass={'build_ext': build_ext},
        ext_modules=[
            Extension('biofrills.cpairutils', ['biofrills/cpairutils.pyx']),
        ],
    )
else:
    setup_args.update(
        ext_modules=[
            Extension('biofrills.cpairutils', ['biofrills/cpairutils.c']),
        ],
    )

setup(**setup_args)

