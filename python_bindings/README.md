# Python Bindings for FBGA (`fbga_py`)

This directory contains the Python bindings for the `FBGA` C++ library. It allows you to use the `FBGA` library directly from Python, providing a convenient interface for working with the library's functionality.

## Installation

To build and install the Python package, you need a C++17 compatible compiler. We provide two recommended ways to set up your environment: using `conda` or using a Python virtual environment (`venv`).

### Option 1: Using Conda (Recommended)

Conda can automatically install a compatible C++ compiler for you, which simplifies the setup process across different platforms (Linux, macOS, Windows).

1. **Create and activate the conda environment:**
    From this directory (`python_bindings/`), run the following command. It will create an environment named `fbga_env` using the provided `environment.yml` file.

    ```bash
    conda env create -f environment.yml
    conda activate vprim_py_env
    ```

2. **Install the package:**
    Once the environment is active, install the package in editable mode using `pip`.

    ```bash
    pip install -e .
    ```

### Option 2: Using `venv` and a system compiler

If you prefer not to use Conda, you can use a standard Python virtual environment. You must ensure you have a C++17 compiler (e.g., GCC, Clang, or MSVC) installed and available in your system's PATH.

```bash
# Navigate to this directory (python_bindings/)
# It is highly recommended to use a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install the package in editable mode (-e).
pip install -e .
```

This command will compile all the C++ source code and create a Python module that you can import.

## Usage

Once installed, you can import and use the classes in your Python code:

```python

```