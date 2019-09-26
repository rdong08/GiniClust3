#!/usr/bin/env python
##################################################
# File Name: setup.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Tue 23 Jul 2019 11:17:18 AM EDT
################################################

#!/usr/bin/env python
from setuptools import setup

setup(
    name='giniclust3',
    version='1.0.0',
    description='Rare cluster identification in single cells',
    url='https://github.com/rdong08/GiniClust3',
    author='Rui Dong',
    author_email='rdong1989@gmail.com',
    install_requires = ['numpy>=1.15.4','scikit-learn>=0.20.2','scipy>=1.2.0','statsmodels>=0.9.0'],
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Software Development :: Build Tools',
        'License :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='GiniClust3 python',
    packages=['giniclust3'],
)
