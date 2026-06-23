# Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study
Code accompanying **Chapter 5** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Paper
Rezaei, M., Fitzgerald, J. G., Wheeler, J. D., Chappell, M. J., and Noel, A. "Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study." To appear in *Proc. IEEE International Conference on Communications (ICC)*, Glasgow, UK, 2026.

*A journal version with extended results is currently in preparation.*

## Overview
This repository implements an analytical diffusion model for molecular signal propagation through multi-layered parallel (Cartesian) channels with heterogeneous transport properties. The framework supports an arbitrary number of layers, flexible transmitter–receiver positioning, and configurable (reflective or permeable) interface boundary conditions, and is applied to a case study modelling molecular transport across the ocean surface microlayer.

## Repository structure
```
src/
├── main.m                      — entry point; sets parameters and runs the simulation
├── greens_function.m           — computes the Green's function for the layered channel
├── boundary_conditions.m       — applies interface/boundary conditions between layers
├── compute_concentration.m     — calculates concentration profile across the parallel channel
└── plot_results.m              — generates output plots

figures/
└── concentration_profile.png   — example output: concentration vs. depth/time

data/
└── ocean_parameters.mat        — diffusion coefficients and layer properties used in the case study
```


## Requirements
- MATLAB R20XX or later
- Toolboxes: [list only if actually required, e.g. none / Symbolic Math Toolbox]

## Usage
1. Open MATLAB and navigate to `src/`
2. Run:
```matlab
   main
```
3. Output: [what it generates — plot? saved file? printed values?]

Key parameters (number of layers, layer thicknesses, diffusion coefficients, transmitter/receiver positions) can be set in `main.m`.

## Example output
[Insert the actual figure — e.g. concentration profile across the ocean surface microlayer from your paper]

## Citation
```bibtex
@inproceedings{rezaei2026multilayer,
  title={Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study},
  author={Rezaei, Mitra and Fitzgerald, J. G. and Wheeler, J. D. and Chappell, M. J. and Noel, A.},
  booktitle={Proc. IEEE International Conference on Communications (ICC)},
  year={2026}
}
```

## Author
Mitra Rezaei — [LinkedIn](https://www.linkedin.com/in/mitra-rezaei-834784159/) · [Thesis hub](https://github.com/Mitra74/phd-thesis-warwick)
