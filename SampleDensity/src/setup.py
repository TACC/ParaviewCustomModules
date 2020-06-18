from setuptools import setup, Extension

setup(ext_modules=[Extension('Sample', ['SampleDensity.c'],),],)
