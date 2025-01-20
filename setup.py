# -*- coding: utf-8 -*-

from os.path import abspath, dirname
from sys import path
from setuptools import setup, find_packages

# To temporarily modify sys.path
SETUP_DIR = abspath(dirname(__file__))

setup(
    name='cobradb',
    version='0.3.0',
    description="""COBRAdb loads genome-scale metabolic models and genome
                   annotations into a relational database.""",
    url='https://github.com/SBRG/cobradb',
    author='Zachary King',
    author_email='zaking@ucsd.edu',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='systems biology, genome-scale model',
    packages=find_packages(),
    install_requires=[
        'SQLAlchemy>=2.0.0,<2.1.0',
        # Be careful upgrading cobra for use with BiGG; ancient models (e.g.
        # iND750.xml) cause errors with cobra versions >=0.15
        'cobra>=0.14.2,<0.30.0',
        'python-libsbml>=5.20.0,<5.21.0',
        'numpy>=2.2.0,<2.3.0',
        'pandas>=2.2.0,<2.3.0',
        'psycopg2>=2.9.0,<2.10.0',
        'biopython>=1.85,<1.86',
        'scipy>=1.15.0,<1.16.0',
        'lxml>=5.3.0,<5.4.0',
        'pytest>=8.3.0,<8.4.0',
        'six>=1.17.0,<1.18.0',
        'escher>=1.8.0,<1.9.0',
        'configparser>=7.1.0,<7.2.0',
    ],
)
