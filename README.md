# Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study
Code accompanying **Chapter 5** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Paper
Rezaei, M., Fitzgerald, J. G., Wheeler, J. D., Chappell, M. J., and Noel, A. "Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study." To appear in *Proc. IEEE International Conference on Communications (ICC)*, Glasgow, UK, 2026.

*A journal version with extended results is currently in preparation.*

## Overview
This repository implements an analytical diffusion model for molecular signal propagation through multi-layered parallel (Cartesian) channels with heterogeneous transport properties. The framework supports an arbitrary number of layers, flexible transmitter–receiver positioning, and configurable (reflective or permeable) interface boundary conditions, and is applied to a case study modelling molecular transport across the ocean surface microlayer.

## Repository structure
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
