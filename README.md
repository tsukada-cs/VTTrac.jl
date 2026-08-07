# VTTrac.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/dev)
[![Build Status](https://github.com/tsukada-cs/VTTrac.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tsukada-cs/VTTrac.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tsukada-cs/VTTrac.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tsukada-cs/VTTrac.jl)

## VTTrac: Velocimetry by Template Tracking
This library provides the julia-language implementation for `VTTrac` (https://github.com/thorinouchi/VTTrac). It does not use module variables, so it should be good for parallel execution.

The algorithm used in this library is the simple template matching of PIV (particle image velocimetry) for monochromatic image-like data, but the matching is conducted multiple times in a Lagrangian manner as in PTV (particle tracking velocimetry) over a number of times specified by the parameter named `ntrac`. The default scoring method for template matching is the cross correlation coefficient, as in the basic PIV. Both forward and backward tracking is available. Use the parameter `itstep`; tracking is backward along time sequence, if it is negative.

## Installation
```julia
using Pkg
Pkg.add("VTTrac")
```

## Quick start

```julia
using VTTrac

# z: image time series, Float32, sized [time, y, x], at least 2 time steps
vtt = VTT(z)

# nsx, nsy: template size; vxhw, vyhw: expected velocity range to search around
setup(vtt, 5, 5; vxhw=2.0, vyhw=2.0, ntrac=5)

# tid0, x0, y0: initial time index / (x, y) position of each template to track
count, status, tid, x, y, vx, vy, score = trac(vtt, tid0, x0, y0)
```

`count[m]` is how many trajectory points succeeded for the `m`-th template, and
`status[m]` says why it stopped otherwise. See the
[How to Use](https://tsukada-cs.github.io/VTTrac.jl/dev/howtouse/) page for a full
walkthrough, including masking and result-screening options.

## Documentation
Documenter.jl generated documentation:
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/dev)

## Contributing
You can contribute to this package by opening issues on GitHub or implementing things yourself and making a pull request. We'd also appreciate more example documents written using `VTTrac`.

## Related packages
* `VTTrac` by Takeshi Horinouchi: https://github.com/thorinouchi/VTTrac
* `pyVTTrac` by Taiga Tsukada: https://github.com/tsukada-cs/pyVTTrac

## References
* Horinouchi, T., S. Tsujino, M. Hayashi, U. Shimada, W. Yanase, A. Wada, and H. Yamada, 2023: Stationary and Transient Asymmetric Features in Tropical Cyclone Eye with Wavenumber-1 Instability: Case Study for Typhoon Haishen (2020) with Atmospheric Motion Vectors from 30-Second Imaging. Monthly Weather Review, 151, 253–273, https://doi.org/10.1175/MWR-D-22-0179.1.
* Tsukada, T., T. Horinouchi, and S. Tsujino, 2024: Wind Distribution in the Eye of Tropical Cyclone Revealed by a Novel Atmospheric Motion Vector Derivation. JGR Atmospheres, 129, e2023JD040585, https://doi.org/10.1029/2023JD040585.
