# Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study

Code accompanying **Chapter 5** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Paper

Rezaei, M., Fitzgerald, J. G., Wheeler, J. D., Chappell, M. J., and Noel, A. "Molecular Communication in Multi-layered Parallel Channels with an Ocean Surface Case Study." To appear in *Proc. IEEE International Conference on Communications (ICC)*, Glasgow, UK, 2026.

*A journal version with extended results is currently in preparation.*

## Overview

This repository implements an analytical diffusion model for molecular signal propagation through multi-layered parallel (Cartesian) channels with heterogeneous transport properties. The mathematical framework presented in the paper supports an arbitrary number of layers, flexible transmitter–receiver positioning, and configurable (reflective or permeable) interface boundary conditions, and is applied to a case study modelling molecular transport across the ocean surface microlayer.

This repository implements a **three-layer** generalization of the paper's two-layer reflective ocean model. In the paper, the sea surface microlayer (SML) is bounded above by a reflective air–water interface and interfaces below with the unbounded bulk ocean. Here, Layer 1 plays the role of the SML — bounded by a reflective (zero-flux) outer boundary, with the point-source transmitter located inside it — Layer 3 plays the role of the unbounded bulk ocean, and Layer 2 is an additional middle layer included to demonstrate the general multi-layer formulation (it is not present in the paper's two-layer reflective case). The coordinate origin and sign of `z` are set up differently from the paper but the underlying physics is the same; see the header of `bounded_three_layer_molecular_diffusion_channel.m` for the explicit geometry mapping.

The model can be extended to an arbitrary number of layers, flexible transmitter positioning (`z0`), and other boundary conditions (e.g. permeable instead of reflective) by following the same pattern used for the existing interfaces — see the extensibility notes in `solveBoundedThreeLayerGreenFunction.m`. The receiver/observation location can be changed by simply changing `z` in the relevant evaluation point, with no further modification needed.

## Repository structure
src/

├── bounded_three_layer_molecular_diffusion_channel.m   — entry point; sets parameters and runs the simulation

├── solveBoundedThreeLayerGreenFunction.m                — assembles and solves the boundary-condition system per frequency

├── evaluateBoundedLayer1GreenFunction.m                 — evaluates the Green's function G in Layer 1 (source/transmitter layer)

├── evaluateBoundedLayer2ImpulseResponse.m               — evaluates the impulse response H2 in Layer 2 (middle layer)

└── evaluateBoundedLayer3ImpulseResponse.m               — evaluates the impulse response H3 in Layer 3 (unbounded outer layer)

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required (uses only core MATLAB functions: `besselj`, `trapz`, `pinv`, `ifft`)

## Usage

1. Open MATLAB and navigate to `src/`
2. Run:
```matlab
   bounded_three_layer_molecular_diffusion_channel
```
3. Output: four figures (combined three-layer response, plus one per layer) showing the time-domain Green's function/impulse response at the specified observation points, console-printed peak values and arrival times, and a `results` struct in the workspace containing the full frequency- and time-domain data.

Key parameters (layer thicknesses, diffusion coefficients, transmitter/observation positions, frequency and Hankel-transform discretization) can be set at the top of `bounded_three_layer_molecular_diffusion_channel.m`.

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
