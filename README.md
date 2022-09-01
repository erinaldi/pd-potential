# Linear confinement in the partially-deconfined phase

Code and data used in the publication _Linear confinement in the partially-deconfined phase_ ([arxiv:2208.14402](https://arxiv.org/abs/2208.14402)) by Gautam, Hanada, Holden, and Rinaldi. 


## Simulation code

The Fortran `F90` code, based on LAPACK, is in the `simulation-code` folder.
The folder includes 2 versions of the code (with and without `OPENMP` directives for multi threading of loops) and a `Makefile`.

## Phase analysis

The python code to do Bayesian analysis of the phase distribution described in Section 4 of the paper is in the `phase-analysis` folder.
The sub-folder `m-values` include the results of the analysis.

## Correlator analysis 

The raw data for the correlator analysis is in the `data-files` folder.
These files are collected from the simulation outputs using the scripts in the `gather-files` folder, and analyzed with the `create-plots` scripts (with output saved in the `products` folder).

## Cite

If you use this code (or parts of it), please consider citing our paper:
```bibtex
@misc{https://doi.org/10.48550/arxiv.2208.14402,
  doi = {10.48550/ARXIV.2208.14402},
  url = {https://arxiv.org/abs/2208.14402},
  author = {Gautam, Vaibhav and Hanada, Masanori and Holden, Jack and Rinaldi, Enrico},
  keywords = {High Energy Physics - Theory (hep-th), High Energy Physics - Lattice (hep-lat), High Energy Physics - Phenomenology (hep-ph), FOS: Physical sciences, FOS: Physical sciences},
  title = {Linear confinement in the partially-deconfined phase},
  publisher = {arXiv},
  year = {2022},
  copyright = {Creative Commons Attribution 4.0 International}
}
```
