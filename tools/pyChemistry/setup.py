################################################################################
# @file setup.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
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
