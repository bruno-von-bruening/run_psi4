# run_psi4

## Install
To run this package please create a python environment through
```bash
conda env create -f install/run_psi4_env.yaml [-n <your-custom-name>]
conda activate run_psi4[|or your custom name]
python -m build  && pip install . # Within repo directory
```

Then run tests:
```bash
pytest # in toplevel or tests/ directory
```
