#!/usr/bin/env python3

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import glob

# Collect all C++ source files from the library's modules (src/utils, src/solvers,
# src/fbga). The paths are relative to this setup.py file.
fbga_sources = (
    glob.glob("../src/utils/*.cc")
    + glob.glob("../src/solvers/*.cc")
    + glob.glob("../src/fbga/*.cc")
)

ext_modules = [
    Pybind11Extension(
        "fbga_py",  # The name of the Python module
        # The source files for the extension: the bindings plus all library sources.
        ["bindings.cc"] + fbga_sources,
        # Include directories: one per module's public include/ dir, matching how
        # each module is consumed via CMake (`#include <utils/...>`, `<solvers/...>`, `<fbga/...>`).
        # pybind11.setup_helpers automatically adds pybind11's includes.
        include_dirs=[
            "../src/utils/include",
            "../src/solvers/include",
            "../src/fbga/include",
        ],
        language='c++',
        cxx_std=20,
    ),
]

setup(
    name="fbga_py",
    version="0.1.0",
    author="Mattia Piazza",
    author_email="mattia.piazza@unitn.it",
    description="Python bindings for the FBGA library",
    long_description="Python bindings for the C++ FBGA library, enabling the use of its functionalities directly from Python.",
    ext_modules=ext_modules,
    # Use the custom build extension command from pybind11
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)