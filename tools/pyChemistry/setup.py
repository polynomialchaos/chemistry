####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
#-----------------------------------------------------------------------------------------------------------------------------------
from setuptools import setup, find_packages

####################################################################################################################################
# Defintions
#-----------------------------------------------------------------------------------------------------------------------------------
version_locals = {}
with open('pychemistry/version.py', 'r') as fp:
    exec(fp.read(), globals(), version_locals)

####################################################################################################################################
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------
setup(
    name='pyChemistry',
    version=str( version_locals['__version__'] ),
    author='Florian Eigentler',
    author_email='f.m.eigentler@gmail.com',
    description='Python package for Chemistry preprocessing',
    packages=find_packages(),
    install_requires=['h5py'],
    entry_points={
        'console_scripts': [
            'pyChemistry=bin.pyChemistry:main',
        ]
    }
)