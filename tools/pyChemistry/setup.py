################################################################################
# @file setup.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from setuptools import setup, find_packages

version_locals = {}
with open('pychemistry/version.py', 'r') as fp:
    exec(fp.read(), globals(), version_locals)

setup(
    name='pyChemistry',
    version=str(version_locals['__version__']),
    author='Florian Eigentler',
    author_email='f.m.eigentler@gmail.com',
    description='Chemistry library preprocessing',
    packages=find_packages(),
    install_requires=['h5py'],
    entry_points={
        'console_scripts': [
            'pyChemistry=pychemistry.bin.pyChemistry:main',
        ]
    }
)
