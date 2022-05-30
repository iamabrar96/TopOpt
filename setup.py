#!/usr/bin/env python
from setuptools import setup

setup(
    name='TopOpt',
    version='0.1.0',
    author='Mohammed Abrar Hyder',
    author_email='ashabrar008@gmail.com',
    description='A Package to implement Topology Optimisation',
    long_description=open('README.md').read(),
    install_requires=["numpy", "torch", "gmsh", "matplotlib", "scipy", "tqdm", "pytest"],

)