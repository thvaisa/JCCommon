import os
import pathlib
from setuptools import setup,find_packages,find_namespace_packages
from Cython.Build import cythonize
from distutils.extension import Extension

setup(
    name = "common",
    version = "0.0.1",
    author = "Timo Väisänen, Aalto university",
    author_email = "",
    description = ("JetCreator common files for Sispo"),
    license = "BSD 2-Clause Simplified License",
    python_requires='>=3.6',
    packages=find_namespace_packages(),
    # Check dependencies
    install_requires=[
        "numpy",
        "orekit",
        "pymesh",
    ],
    include_package_data=True,
)
