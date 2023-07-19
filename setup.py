#!/usr/bin/env python

from setuptools import setup, find_packages

# Read requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='oceantracker',
      python_requires='>=3.10,<3.11',
      version='0.4.0',
      description='Fast offline Lagrangian particle tracking in the Ocean',
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      author='Ross Vennell',
      author_email='ross.vennell@cawthron.org.nz',
      url='https://oceantracker.github.io/oceantracker/',
      packages=find_packages(),
      install_requires=requirements,
)
