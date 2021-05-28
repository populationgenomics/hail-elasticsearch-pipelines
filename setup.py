#!/usr/bin/env python

"""
Setup script for the Python package
- Used for development setup with `pip install --editable .`
- Parsed by conda-build to extract version and metainfo
"""

import setuptools
import os
from os.path import join, relpath

PKG = 'seqr-loader'


def find_package_files(dirpath, package, skip_exts=None):
    paths = []
    for (path, _dirs, fnames) in os.walk(join(package, dirpath)):
        for fname in fnames:
            if skip_exts and any(fname.endswith(ext) for ext in skip_exts):
                continue
            fpath = join(path, fname)
            paths.append(relpath(fpath, package))
    return paths


setuptools.setup(
    name='seqr-loader',
    # This tag is automatically updated by bump2version
    version='1.0.5',
    description='The hail scripts in this repo can be used to pre-process variant callsets and export them to elasticsearch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/hail-elasticsearch-pipelines',
    license='MIT',
    packages=['hail_scripts', 'model'],
    package_dir={'model': 'batch_seqr_loader/model', 'hail_scripts': 'hail_scripts'},
    include_package_data=True,
    zip_safe=False,
    scripts=['batch_seqr_loader/seqr_load.py'],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
