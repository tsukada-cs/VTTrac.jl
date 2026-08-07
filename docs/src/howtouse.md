```@meta
CurrentModule = VTTrac
```

# How to Use

## Overview

Using VTTrac generally follows five steps:

1. Wrap your image time series in a [`VTT`](@ref) object.
2. Configure tracking parameters with [`setup`](@ref).
3. Choose the initial time/position of each template you want to track.
4. Run [`trac`](@ref).
5. Interpret the returned trajectories, using `count`/`status` to tell which points succeeded.

## 1. Prepare the data

`z` must be a 3-D `Float32` array in `[time, y, x]` order, with at least 2 time steps.
`t` is optional and defaults to `0:size(z,1)-1`; it only needs to be supplied if your
time steps are irregularly spaced or you want physically-dimensional velocities.

```julia
using VTTrac

nt, ny, nx = 10, 100, 100
tax = Vector{Float64}([0:nt-1;])
xax = [0:nx-1;]
yax = [0:ny-1;]
xg = permutedims(repeat(xax', outer=(length(yax),1,length(tax))), (3,1,2))
yg = permutedims(repeat(yax, outer=(1,length(xax),length(tax))), (3,1,2))
tg = repeat(tax, outer=(1, length(yax), length(xax)))
k = 2pi/10
cx, cy = 1.2, 1.2
z = sin.(k*(xg .- cx*tg)) .* cos.(k*(yg .- cy*tg))
z = Array{Float32,3}(z) # a wave pattern drifting at (vx, vy) = (cx, cy)
```

## 2. Create the `VTT` object

```julia
vtt = VTT(z, t=tax)
```

If some pixels should never be used for scoring (e.g. sensor gaps, quality-flagged
pixels), pass a `Bool` array of the same size as `z` via `mask=...`, where `true` means
"ignore this pixel". See [Masking invalid pixels](@ref) below.

## 3. Configure tracking with `setup`

```julia
setup(vtt, 5, 5; vxhw=1.8, vyhw=1.8, ntrac=5, subgrid=true, score_method="xcor")
```

- `nsx`, `nsy` (positional): the template sub-image size.
- Either `vxhw`/`vyhw` (an expected velocity range) **or** `ixhw`/`iyhw` (a pixel search
  half-width) is mandatory — they determine how far around the first guess (or the
  previous step's result) the next match is searched for.
- `ntrac`: how many times each initial template is tracked forward (or backward, see
  `itstep` below). Use `ntrac >= 2` if you want velocity-change screening (`vxch`/`vych`).
- `score_method`: `"xcor"` (cross-correlation, the default) or `"ncov"` (normalized
  covariance) — see [Available scoring methods](@ref) on the Home page.
- `itstep`: the index step between tracked frames; negative values track backward in
  time.
- `subgrid`, `subgrid_gaus`: whether/how to refine the matched position to sub-pixel
  accuracy (elliptic-paraboloid by default, Gaussian if `subgrid_gaus=true`).

See [`setup`](@ref) for the full parameter list, including the result-screening options
covered in [Screening options](@ref) below.

## 4. Choose initial points and track

`tid0`, `x0`, `y0` locate each template's starting time index and (x, y) position
(index-based; non-integer positions are allowed and are read out via bilinear
interpolation). They must share the same shape — any shape works, not just vectors.

```julia
n = 6
tid0 = ones(Int, n)
x0 = 1 .+ [0:n-1;]*2.5 .+ 7.5
y0 = 1 .+ [0:n-1;]*1.0 .+ 10.5

count, status, tid, x, y, vx, vy, score, zss, score_ary = trac(vtt, tid0, x0, y0)
```

`out_subimage=true`/`out_score_ary=true` additionally return the template sub-images and
full sliding-score arrays along the trajectory, which is useful for diagnosing why a
particular point stopped tracking (see `zss`, `score_ary` in [`trac`](@ref)).

## 5. Interpret the results

Every output's first dimension (or, for `count`/`status`, the whole array) is shaped
like `tid0`/`x0`/`y0`. `tid`, `x`, `y` have `ntrac+1` rows — the initial point plus one
row per successful step; `vx`, `vy`, `score` have `ntrac` rows, one per step.

- `count[m]`: how many trajectory points (including the initial one) are valid for
  the `m`-th template — `ntrac+1` means every step succeeded.
- `status[m]`: *why* tracking stopped for the `m`-th template (`0` if it never
  stopped early). See the table below.
- By default (`to_missing=true`), positions past `count[m]` are filled with
  `missing` rather than the raw sentinel values (`vtt.fmiss`/`vtt.imiss`).

| `status` | Meaning |
|:--------:|---------|
| `0` | Completed successfully (or the step hasn't been reached yet) |
| `1` / `5` | The start/end time index for this step was out of range |
| `2` | The template sub-image could not be read (out of the image bounds, or missing/masked data) |
| `3` | The template sub-image failed the `Cth` contrast check |
| `4` | The template sub-image failed the `peak_inside_th` check |
| `6` | The sliding score could not be computed (out of bounds, or missing/masked data) |
| `7` | No valid score peak was found, or the peak fell on the edge of the search window |
| `8` | The peak score was below `Sth0` (1st step) or `Sth1` (subsequent steps) |
| `9` | The velocity change along the trajectory exceeded `vxch`/`vych` |

## Masking invalid pixels

Pass `mask` (same shape as `z`, `true` = ignore) when constructing `VTT` to skip
untrustworthy pixels during scoring:

```julia
mask = falses(size(z))
mask[:, 1:5, :] .= true # e.g., the top 5 rows are known to be bad

vtt = VTT(z, t=tax, mask=mask)
setup(vtt, 5, 5; vxhw=1.8, vyhw=1.8, ntrac=5, min_samples=10)
```

`min_samples` sets how many jointly-unmasked pixels are required (between the template
and the candidate window) before a score is computed there at all; positions with fewer
are treated as undefined, the same as `status = 6`.

## Screening options

- `Sth0`, `Sth1` (defaults `0.8`, `0.7`): minimum score required for the 1st and
  subsequent steps, respectively.
- `vxch`, `vych`: if **both** are set (positive), a step is rejected when the velocity
  changes by more than this amount from the previous step (`status = 9`); this needs
  `ntrac >= 2` to have any effect.
- `peak_inside_th`: if set, an initial template is only used when it has a prominent
  interior peak or trough (rejects flat/edge-dominated templates, `status = 4`).
- `Cth`: if set, an initial template is only used when its max−min exceeds this value
  (rejects low-contrast templates, `status = 3`).
- `use_init_temp`: if `true`, the template read out at the very first step is reused for
  every subsequent step, instead of re-reading a fresh template around the current
  position at each step.

## Tracking backward in time

Set `itstep` to a negative value in [`setup`](@ref) to track backward along the time
axis from `tid0`.
