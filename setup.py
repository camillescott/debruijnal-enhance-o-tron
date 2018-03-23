#!/usr/bin/env python
import os

from setuptools import setup
from setuptools import find_packages

__path__ = os.path.dirname(__file__)
__version__ = open(os.path.join(__path__, 'debruijnal_enhance_o_tron', 'VERSION')).read().strip()

packages = find_packages(exclude=['docs'])

setup(
    name = 'debruijnal-enhance-o-tron',
    version = __version__,

    author = 'Camille Scott',
    author_email = 'camille.scott.w@gmail.com',

    url = 'https://github.com/camillescott/debruijnal-enhance-o-tron',
    description = 'pytest fixtures for testing de Bruin Graphs',
    license = 'BSD-3-Clause',

    packages = packages,
    zip_safe = False,
    platforms = 'any',

    install_requires = ['pytest'],
    tests_require = ['pytest', 'pytest-runner'],
    setup_requires = ['pytest', 'pytest-runner'],

    keywords = 'pytest bioinformatics testing',
    classifiers = [
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Testing',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
