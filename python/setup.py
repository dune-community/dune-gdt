#!/usr/bin/env python

import sys
from setuptools import setup

setup(name='dune.gdt',
      version='0.3-dev',
      namespace_packages=['dune'],
      description='Python for Dune-Gdt',
      author='The dune-gdt devs',
      author_email='dune-gdt-dev@listserv.uni-muenster.de',
      url='https://github.com/dune-community/dune-gdt',
      packages=['dune.gdt'],
      install_requires=['jinja2', 'where'],
      )
