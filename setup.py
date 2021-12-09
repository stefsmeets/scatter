#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path

# www.pythonhosted.org/setuptools/setuptools.html

setup(
    name="scatter",
    version="1.0.2",
    description="A python tool to plot and output atomic scattering factors",

    author="Stef Smeets",
    author_email="s.smeets@esciencecenter.nl",
    url="https://github.com/stefsmeets/scatter",

    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'License :: OSI Approved :: MIT License',
    ],

    packages=["scatter"],

    install_requires=["numpy", "matplotlib"],

    package_data={
        "": ["LICENCE",  "readme.md", "setup.py", "scatter.py"],
    },

    entry_points={
        'console_scripts': [
            'scatter = scatter.scatter:main',
        ]
    }

)
