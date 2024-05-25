# VTTrac.jl 

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tsukada-cs.github.io/VTTrac.jl/dev)
[![Build Status](https://github.com/tsukada-cs/VTTrac.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tsukada-cs/VTTrac.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tsukada-cs/VTTrac.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tsukada-cs/VTTrac.jl)

Documentation for [VTTrac](https://github.com/tsukada-cs/VTTrac.jl).

## VTTrac: Velocimetry by Template Tracking
This library provides the julia-language implementation for `VTTrac` (https://github.com/thorinouchi/VTTrac). It does not use module variables, so it should be good for parallel execution.

The algorithm used in this library is the simple template matching of PIV (particle image velocimetry) for monochromatic image-like data, but the matching is conducted multiple times in a Lagrangian manner as in PTV (particle tracking velocimetry) over a number of times specified by the parameter named `ntrac`. The default scoring method for template matching is the cross correlation coefficient, as in the basic PIV. Both forward and backward tracking is available. Use the parameter `itstep`; tracking is backward along time sequence, if it is negative.