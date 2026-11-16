# run_psi4

This package builds input files for <tt>psi4</tt>[^1].
Automatic basis set generation, for e.g. doubly-augmented functions not included in <tt>psi4</tt> are implemented
through the <tt>basis-set-exchange</tt> package[^2].

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

# References
[^1]: Smith, D. G. A.; Burns, L. A.; Simmonett, A. C.; Parrish, R. M.; Schieber, M. C.; Galvelis, R.; Kraus, P.; Kruse, H.; Di Remigio, R.; Alenaizan, A.; James, A. M.; Lehtola, S.; Misiewicz, J. P.; Scheurer, M.; Shaw, R. A.; Schriber, J. B.; Xie, Y.; Glick, Z. L.; Sirianni, D. A.; O’Brien, J. S.; Waldrop, J. M.; Kumar, A.; Hohenstein, E. G.; Pritchard, B. P.; Brooks, B. R.; Schaefer, H. F.; Sokolov, A. Yu.; Patkowski, K.; DePrince, A. E.; Bozkaya, U.; King, R. A.; Evangelista, F. A.; Turney, J. M.; Crawford, T. D.; Sherrill, C. D. P SI4 1.4: Open-source software for high-throughput quantum chemistry. *The Journal of Chemical Physics* **2020**, 152 (18), 184108. https://doi.org/10.1063/5.0006002.
See github project: https://github.com/psi4/psi4/
[^2]: Pritchard, B. P.; Altarawy, D.; Didier, B.; Gibson, T. D.; Windus, T. L. New Basis Set Exchange: an open, Up-to-Date resource for the molecular sciences community. *Journal of Chemical Information and Modeling* **2019**, 59 (11), 4814–4820. https://doi.org/10.1021/acs.jcim.9b00725.
See github project: https://github.com/MolSSI-BSE/basis_set_exchange
