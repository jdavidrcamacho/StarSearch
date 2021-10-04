#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup

setup(name='starsearch',
      version='0.3',
      description='Package to dig into the ESO archives',
      author='João Camacho',
      author_email='joao.camacho@astro.up.pt',
      license='MIT',
      url='https://github.com/jdavidrcamacho/starsearch',
      packages=['starsearch'],
      install_requires=[
        'numpy',
        'astroquery',
        "astropy",
      ],
     )
