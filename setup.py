from setuptools import setup
import os
import sys
from os import path


__version__ = '0.0.1'


def readme():
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
    return long_description


setup(
    name = 'pygff',
    version = __version__,
    description = 'Utility program for parsing GFF3 files',
    long_description = readme(),
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/betteridiot/pygff',
    author = 'Marcus D. Sherman',
    author_email = 'mdsherman@betteridiot.tech',
    license = 'BSD 3-Clause',
    install_requires=[
        'numpy',
        'pandas',
    ],
    packages = ['pygff', 'pygff.tests'],
    package_data = {'pygff': ['LICENSE', 'CODE_OF_CONDUCT.md', 'pygff/tests/data/*']}
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis'
    ],
    keywords = 'gff parsing bioinformatics genomics',
    include_package_data = True,
    zip_safe = False
)   
     
