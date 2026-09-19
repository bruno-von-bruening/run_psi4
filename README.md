# run_psi4

## 0.01.00 Start of package
Managing psi4 as a package (transfered form previous script based implementation)

## Install
To run this package please create a python environment through
```python
conda env create -f install/run_psi4_env.yaml [-n <your-custom-name>]
conda activate run_psi4[|or your custom name]
python -m build  && pip install . # Within repo directory
```

## Current branch
Make new parser and include settings for psi4

Achieve that through reading in a config model

Note that linear response does not work for mp2 and ccsd(t) (but what is here computed anyway?)
