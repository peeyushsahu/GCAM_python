#!/usr/bin/env python

"""Description

Setup script for GCAM -- Gene Celltype Association Miner

Copyright (c) 2015 Peeyush Sahu <peeyush.sahu@imt.uni-marburg.de>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  beta
@version: $Revision$
@author:  Peeyush Sahu
@contact: peeyush.sahu@imt.uni-marburg.de
"""


import sys
from setuptools import setup, Extension


def main():
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    setup(name="GCAM",
          version="1.0.0",
          description="Gene Celltype Association Miner",
          author='Peeyush Sahu',
          author_email='peeyush.sahu@imt.uni-marburg.de',
          url='https://github.com/peeyushsahu/GCAM_python.git',
          packages=['GCAM'],
          package_dir={'GCAM' : 'GCAM'},
          package_data={'GCAM': ['resources/*']},
          scripts=['bin/GCAM'],

          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: Linux',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python',
          ],
          install_requires=[
              'numpy>=1.6',
              'matplotlib',
              'pandas',
              #'biopython',
          ],
          )

if __name__ == '__main__':
    main()

