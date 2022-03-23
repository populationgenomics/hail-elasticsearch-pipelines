#!/usr/bin/env python

"""
Setup script for the Python package
- Used for development setup with `pip install --editable .`
- Parsed by conda-build to extract version and metainfo
"""

import setuptools

PKG = 'seqr-loader'

setuptools.setup(
    name='seqr-loader',
    # This tag is automatically updated by bump2version
    version='1.2.0',
    description='The hail scripts in this repo can be used to pre-process variant callsets and export them to elasticsearch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/hail-elasticsearch-pipelines',
    license='MIT',
    packages=['hail_scripts', 'lib', 'lib.model'],
    package_dir={
        'lib': 'luigi_pipeline/lib',
        'lib.model': 'luigi_pipeline/lib/model',
    },
    include_package_data=True,
    zip_safe=False,
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
