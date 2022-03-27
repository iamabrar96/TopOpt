# Installation of the Program

This document details the manual installation of the Program.

### Python version

The program itself is not bound to a certain python version, as long as the package requirements for the program can be met by the
version of your choice. The program has been successfully tested for python 3.8.5.

### Python packages

In order to run the program you have to have the following packages installed:

* torch (PyTorch)
* numpy
* scipy
* gmsh
* matplotlib
* tqdm

You can install each Python package with

```sh
$ pip install libraryname
```

or just install all the libraries

```sh
$ pip install -e .
```