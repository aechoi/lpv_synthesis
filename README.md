# lpv_synthesis

An implementation of robust control synthesis for LPV systems with IQC defined uncertainties. Details can be found in the paper "Robust Synthesis for Linear Parameter Varying Systems Using Integral Quadratic Constraints" by S. Wang, H. Pfifer, and P. Seiler.

# dev guide

Do the following commands in this directory.

- Upgrade pip `python -m pip install --upgrade pip`
- Initialize a python virtual environment `python -m venv [virtual environment name]`
- Activate the virtual environment
  - Windows `.\[virtual env name]\Scripts\activate`
  - Linux\OSX `source .\[virtual env name]\bin\activate`
- Install requirements `python -m pip install -r requirements.txt`
- Install this package in editable mode (for pytest and relative imports) `pip install -e .`

The structure of `src\lpv_synthesis` should be identical to the structure of `tests\`, following the convention found [here](https://docs.pytest.org/en/stable/explanation/goodpractices.html#tests-outside-application-code).

