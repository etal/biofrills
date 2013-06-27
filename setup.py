#!/usr/bin/env python

"""Bioinformatics utilities for molecular sequence analysis."""

from glob import glob
from os.path import dirname, join

DIR = (dirname(__file__) or '.')

setup_args = dict(
    name='biofrills',
    version='0.3',
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@gmail.com',
    url='http://github.com/etal/biofrills',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['biofrills', 'biofrills.stats'],
    #scripts=glob(join(DIR, 'scripts', '*'))
)

try:
    from setuptools import setup
    from setuptools.extension import Extension
    # Dependencies for easy_install:
    setup_args.update(
        install_requires=[
            'biopython >= 1.60',
            'Cython >= 0.15',
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

