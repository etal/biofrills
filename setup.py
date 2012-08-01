#!/usr/bin/env python

"""Bioinformatics utilities for sequence analysis."""

setup_args = {}

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

setup_args.update(
    name='biofrills',
    version='0.0.0-dev',
    description=__doc__,
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://etalog.blogspot.com/',
    packages=['biofrills',
             ],
)

try:
    from os.path import dirname
    from Cython.Distutils import build_ext
    setup_args.update(
        cmdclass={'build_ext': build_ext},
        ext_modules=[
            Extension('biofrills.cpairutils',
                      [(dirname(__file__) or '.') +
                       '/biofrills/cpairutils.pyx']),
        ],
    )
except ImportError:
    pass

setup(**setup_args)

