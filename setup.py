#!/usr/bin/env python

"""Bioinformatics utilities for sequence analysis."""

from glob import glob
from os.path import dirname, join

DIR = (dirname(__file__) or '.')

setup_args = dict(
    name='biofrills',
    version='0.0.0',
    description=__doc__,
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://github.com/etal/biofrills',
    packages=['biofrills'],
    #scripts=glob(join(DIR, 'scripts', '*'))
)

try:
    from setuptools import setup
    from setuptools.extension import Extension
    # Dependencies for easy_install:
    setup_args.update(
        install_requires=[
            'Biopython >= 1.59',
            'Cython >= 0.11',
        ])
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    setup_args.update(
        cmdclass={'build_ext': build_ext},
        ext_modules=[
            Extension('biofrills.cpairutils',
                      [join(DIR, 'biofrills', 'cpairutils.pyx')]),
        ],
    )
except ImportError:
    pass

setup(**setup_args)

