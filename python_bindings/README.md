# Python Bindings for FBGA (`fbga_py`)

This directory contains the Python bindings for the `FBGA` C++ library. It allows you to use the `FBGA` library directly from Python, providing a convenient interface for working with the library's functionality.

## Installation

To build and install the Python package, you need a C++20 compatible compiler. We provide two recommended ways to set up your environment: using `conda` or using a Python virtual environment (`venv`).

### Option 1: Using Conda (Recommended)

Conda can automatically install a compatible C++ compiler for you, which simplifies the setup process across different platforms (Linux, macOS, Windows).

1. **Create and activate the conda environment:**
    From this directory (`python_bindings/`), run the following command. It will create an environment named `fbga_py_env` using the provided `fbga_env.yaml` file.

    ```bash
    conda env create -f fbga_env.yaml
    conda activate fbga_py_env
    ```

2. **Install the package:**
    Once the environment is active, install the package in editable mode using `pip`.

    ```bash
    pip install -e .
    ```

### Option 2: Using `venv` and a system compiler

If you prefer not to use Conda, you can use a standard Python virtual environment. You must ensure you have a C++20 compiler (e.g., GCC, Clang, or MSVC) installed and available in your system's PATH.

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
import fbga_py as fb

def gg_upper(ay, v):
    return 10.0  # replace with a real friction-circle upper bound

def gg_lower(ay, v):
    return -10.0

gg_range = fb.GgRangeMaxMin()
gg_range.min = lambda v: -10.0
gg_range.max = lambda v: 10.0

solver = fb.Fb2d(gg_upper, gg_lower, gg_range)
total_time = solver.compute(SS=[0.0, 10.0, 20.0], KK=[0.0, 0.0, 0.0], v0=20.0)
```

See `examples/example_fb2d_plots.py` for a complete, runnable example with plotting.

### fbga3d (`FbgaIndy` / `FbgaMoto`)

The 3D forward-backward solver is also bound, generic over the vehicle's GG-diagram model
(`GggvIndy`, spline-based; `GggvMoto`, closed-form). Both `FbgaIndy` and `FbgaMoto` expose the
same API (`compute`, `get_solution`, `eval_V`, `eval_A_tilde_x`, ...) since they're both
instantiations of the same C++ template -- see `src/fbga3d/include/fbga3d/fbga3d_solver.hh`.

```python
import fbga_py as fb

TOA = fb.TrajectoryOffsetAndAnglesContainer(
    offset=fb.TrajectoryOffsetContainer(n=[0.0, 0.0, 0.0], chi=[0.0, 0.0, 0.0]),
    reference=fb.RoadAnglesAndDerivativesContainer(
        mu=[0.0, 0.0, 0.0], phi=[0.0, 0.0, 0.0], theta=[0.0, 0.05, 0.1],
        mu_prime=[0.0, 0.0, 0.0], phi_prime=[0.0, 0.0, 0.0], theta_prime=[0.001, 0.005, 0.0],
        abscissa=[0.0, 50.0, 100.0],
    ),
    adherence=fb.AdherenceContainer(alpha=[1.0, 1.0, 1.0]),
)

solver = fb.FbgaMoto()  # or fb.FbgaIndy() -- reads ./data/INDY/*.npy relative to the CWD
total_time = solver.compute(TOA, v_initial=20.0)
sol = solver.get_solution()
```

Plotting (`eval_shell_plot*`) and constraint-satisfaction-checking helpers from the original
FBGA_3D aren't ported to the C++ core yet (see `.claude/FBGA3D_INTEGRATION_PLAN.md`), so they
aren't bound here either.
