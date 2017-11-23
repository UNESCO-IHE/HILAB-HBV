# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 08:59:02 2017

@author: chaco3
"""

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sourcefiles = ['HBV96x.pyx', ]

extensions = [Extension("HBV96x", sourcefiles)]

setup(
    ext_modules = cythonize(extensions)
)