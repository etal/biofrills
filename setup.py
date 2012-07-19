#!/usr/bin/env python

try:
    from setuptools import setup
    from setuptools.extension import Extension
except ImportEror:
    from distutils.core import setup
    from distutils.extension import Extension

setup_args = dict(
    name='biofrills',
    version='0.0.0-dev',
    description='Bio scripts, parsers and utilities',
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://etalog.blogspot.com/',
    packages=['biofrills',
             ],
)

try:
    from Cython.Distutils import build_ext
    setup_args.update(
        cmdclass={'build_ext': build_ext},
        ext_modules=[
            Extension('biofrills.cpairutils', ['biofrills/cpairutils.pyx']),
        ],
    )
except ImportError:
    pass

setup(**setup_args)

