#!/usr/bin/env python

from setuptools import setup

setup(
    name='oasa',
    version='0.14.0',
    description="OASA is a free cheminformatics library written in Python",
    author="Beda Kosata",
    author_email="beda@zirael.org",
    maintainer="Adel Qalieh",
    maintainer_email="adelq@med.umich.edu",
    url="http://bkchem.zirael.org/oasa_en.html",
    long_description=open("README.rst").read(),
    packages=['oasa', 'oasa/graph'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization"
    ]
)
