module VTTrac
export VTT, setup, trac

using Printf
using Statistics

import Logging

mutable struct VTT
    # data on which tracking is made
    nx::Int # image size x
    ny::Int # image size y
    nt::Int # time length
    z::Array{Float32,3}
    visible::Union{BitArray{3}, Array{Bool,3}}
    t::Vector{Float64}
    dtmean::Float64
    zmiss::Float32
    fmiss::Float64
    imiss::Int

    # tracking parameters
    chk_zmiss::Bool # if true, check missing values in `z` (image)
    chk_mask::Bool # if true, check mask in `mask` (mask of image)
    nsx::Int # sub-image size x
    nsy::Int # sub-image size y
    vxhw::Float64  # velocities corresponding to `ixhw` through `dtmean`
    vyhw::Float64  # velocities corresponding to `iyhw` through `dtmean`
    ixhw::Int # max displacement x for template matching
    iyhw::Int # max displacement y for template matching
    vxch::Float64
    vych::Float64
    itstep::Int
    ntrac::Int

    subgrid::Bool
    subgrid_gaus::Bool #true: subgrid peak finding is by gaussian; false: e-paraboloid
    score_method::String
    score_th0::AbstractFloat
    score_th1::AbstractFloat
    chk_peak_inside::Bool
    peak_inside_th::Float32 #threshold for the peak-inside screening(unused if<0)
    chk_min_contrast::Bool
    min_contrast::Float32 # minimum contrast in the template (unused if=<0)
    use_init_temp::Bool # if true, always use initial template submimages
    setuped::Bool
    min_visible::Int # minimum number of visible values to calculate score when `chk_mask` is true.

    """
        VTT(z[, t, visible, zmiss, fmiss, imiss])
    
    Sets data for tracking (you need to set parameters separately).
    
    # Arguments
    - `z::Array{Float32,3}`: Array of image-like data (in dimensions [time, y, x]). `z[i]` contains `i`-th image data.
    - `t::Vector{Float64}`: Times at which the images are for.
    - `zmiss::Union{Real, Nothing}=nothing`: Missing value used in `z`.
    - `mask::Union{Array{Bool,3}, Nothing}=nothing`: Mask to ignore when calculate score (true positions are ignored).
    - `fmiss::Real=-999.0`: Missing value to be set for Real.
    - `imiss::Integer=-999`: Missing value to be set for Integer.
    """
    function VTT(z::Array{Float32,3}; t::Union{Vector{Float64}, Nothing}=nothing, mask::Union{BitArray{3}, Array{Bool,3}, Nothing}=nothing, zmiss::Union{Real, Nothing}=nothing, fmiss::Real=-999.0, imiss::Int=-999)
        o = new()
        o.z = z
        o.nt, o.ny, o.nx = size(z)
        if isnothing(t)
            t = Vector{Float64}([0:o.nt-1;])
        end
        size(t) !== (o.nt,) && throw(ArgumentError("`size(t)` must be `(size(z)[begin],1)`"))
        o.t = t
        o.dtmean = (t[end]-t[begin])/(o.nt-1)

        if isnothing(zmiss)
            o.chk_zmiss = false
            o.zmiss = Float32(0.0)
        else
            o.chk_zmiss = true
            o.zmiss = Float32(zmiss)
        end
        o.fmiss = fmiss
        o.imiss = imiss

        if isnothing(mask) || !any(mask)
            o.chk_mask = false
        else
            size(mask) != size(z) && throw(ArgumentError("`size(mask)` must be `(size(z)`"))
            o.chk_mask = true
            o.visible = .!mask
        end

        o.setuped = false
        return o
    end
end

"""
    setup(o, nsx, nsy, vxhw, vyhw, [ixhw, iyhw, subgrid, subgrid_gaus,
        itstep, ntrac, score_method, score_th0, score_th1, vxch, vych,
        peak_inside, peak_inside_th, min_contrast, use_init_temp, min_visible])

Setup for tracking.

# Arguments
- `o::VTT`: The object.
- `nsx::Integer, nsy::Integer`: Submimage x & y sizes (x:1st, y:2nd dim).
- `vxch::Union{Real, Nothing}, vyhw::Union{Real, Nothing}`: (either `v[xy]hw` or `i[xy]hw` are MANDATORY).
    the dimensions along which to perform the computation.
    search velocity range half sizes to set `i[xy]hw`.
    Seach at least to cover +-v?hw around the first guess or previous step.
    (the result can be outside the range.)
- `ixhw::Union{Int, Nothing}, iyhw::Union{Int, Nothing}`: (either `v[xy]hw` or `i[xy]hw` are MANDATORY)
    Max displacement fro template match (can be set indirecly through `v[xy]hw`).
- `subgrid::Bool=true`: Whether to conduct subgrid tracking.
- `subgrid_gaus::Bool=true`: Whether subgrid peak finding is by gaussian.
- `itstep::Integer=1`: Step of `t`'s used (skip if >1).
- `ntrack::Integer=2`: Max tracking times from initial loc.
- `score_method::String="xcor"`: `"xcor"` for cross-correlation, `"ncov"` for normalized covariance.
- `score_th0::AbstractFloat=0.8`: The minimum score required for the 1st tracking.
- `score_th1::AbstractFloat=0.7`: The minimum score required for subsequent tracking.
- `vxch::Union{Real, Nothing}=nothing`: If non-`nothing`, the max tolerant vx
    change between two consecutive tracking.
- `vych::Union{Real, Nothing}=nothing`: If non-`nothing`, the max tolerant vy
    change between two consecutive tracking.
- `peak_inside_th::Union{Real, Nothing}=nothing`: If non-`nothing`, an initial template is used only when
    it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.
- `min_contrast::Union{Real, Nothing}=nothing`: If non-`nothing`, an initial template is used only when 
    it has a difference in max and min greater than its value.
- `min_visible::Int=1`: Minimum number of visible values to calculate score when `chk_mask` is true.
"""
function setup(o::VTT, nsx::Int, nsy::Int; vxhw::Union{Real, Nothing}=nothing, vyhw::Union{Real, Nothing}=nothing,
            ixhw::Union{Int, Nothing}=nothing, iyhw::Union{Int, Nothing}=nothing, subgrid::Bool=true,
            subgrid_gaus::Bool=false, itstep::Int=1, ntrac::Int=2, score_method::String="xcor",
            score_th0::AbstractFloat=0.8, score_th1::AbstractFloat=0.7, vxch::Union{Real, Nothing}=nothing,
            vych::Union{Real, Nothing}=nothing, peak_inside_th::Union{Real, Nothing}=nothing,
            min_contrast::Union{Real, Nothing}=nothing, use_init_temp::Bool=false, min_visible::Int=1)
    o.nsx = nsx
    o.nsy = nsy
    if vxhw !== nothing
        ixhw !== nothing && throw(ArgumentError("`v[xy]hw` and `i[xy]hw` must not be set simultaneously"))
        vyhw === nothing && throw(ArgumentError("vxhw and vyhw must be set simultaneously"))
        set_ixyhw_from_v!(o, vxhw, vyhw)
    elseif ixhw !== nothing
        iyhw === nothing && throw(ArgumentError("ixhw and iyhw must be set simultaneously"))
        set_ixyhw_directly!(o, ixhw, iyhw)
    else
        throw(ArgumentError("either `i[xy]hw` or `v[xy]hw` must be specified"))
    end

    vxch === nothing && (vxch = -999.0) # <=0 for nothing (not to set)
    vych === nothing && (vych = -999.0) # <=0 for nothing (not to set)

    if isnothing(peak_inside_th)
        peak_inside_th = -1.0  # negative, meaning unused
    end
    peak_inside_th = Float32(peak_inside_th)
    if isnothing(min_contrast)
        min_contrast = -1.0   # negative, meaning unused
    end
    min_contrast = Float32(min_contrast)
    
    set_basic!(o, nsx, nsy, itstep, ntrac)
    set_optional!(o, subgrid, subgrid_gaus, score_method, score_th0, score_th1, peak_inside_th, min_contrast, vxch, vych, use_init_temp, min_visible)
    o.setuped = true
end

function set_nsx(o::VTT, value::Int)
    o.nsx = value
end
function set_nsy(o::VTT, value::Int)
    o.nsy = value
end
function set_itstep(o::VTT, value::Int)
    o.itstep = value
end
function set_ntrac(o::VTT, value::Int)
    o.ntrac = value
end
function set_subgrid(o::VTT, value::Bool)
    o.subgrid = value
end
function set_subgrid_gaus(o::VTT, value::Bool)
    o.subgrid_gaus = value
end
function set_score_method(o::VTT, value::String)
    o.score_method = value
end
function set_score_th0(o::VTT, value::Real)
    o.score_th0 = value
end
function set_score_th1(o::VTT, value::Real)
    o.score_th1 = value
end
function set_peak_inside_th(o::VTT, value::Real)
    o.peak_inside_th = value
end
function set_min_contrast(o::VTT, value::Real)
    o.min_contrast = value
end
function set_use_init_temp(o::VTT, value::Bool)
    o.use_init_temp = value
end

"""
    set_ixyhw_from_v!(o, vxch, vyxh)

Sets the tracking parameters `i[xy]hw` from velocities (v[xy]hh).

# Arguments
- `o::VTT`: The object.
- `vxhw::Float64`: The range over which vx is searched around initial guess.
- `vyhw`::Float64: The range over which vy is searched around initial guess.
```

"""
function set_ixyhw_from_v!(o::VTT, vxhw::Float64, vyhw::Float64)
    o.vxhw = vxhw
    o.vyhw = vyhw
    o.ixhw = ceil(abs(vxhw * o.dtmean)) + 1 # max displacement
    o.iyhw = ceil(abs(vyhw * o.dtmean)) + 1 # +1 is margin to find peak
end

"""
    set_ixyhw_from_v!(o, ixch, iyxh)

Sets the tracking parameters `i[xy]hw`.

# Arguments
- `o::VTT`: The object.
- `ixhw::Float64`: The range over which next x is searched around initial guess.
- `iyhw`::Float64: The range over which next y is searched around initial guess.
"""
function set_ixyhw_directly!(o::VTT, ixhw::Int, iyhw::Int)
    o.ixhw = ixhw
    o.iyhw = iyhw
    o.vxhw = ixhw/o.dtmean - 1 # max displacement
    o.vyhw = iyhw/o.dtmean - 1 # -1 is from margin to find peak
end

"""
    set_basic!(o, nsx, nsy, itstep, ntrac)

Sets basic (mandatory) tracking parameters. Also sets default vals for optional params.

# Arguments
- `o::VTT`: The object.
- `nsx::Integer`: The template subimage size (x).
- `nsy::Integer`: The template subimage size (y).
- `itstep::Integer`: Index-based time step for tracking (can be negative).
- `ntrac::Integer`: Number of times for each initial template is tracked.
"""
function set_basic!(o::VTT, nsx::Int, nsy::Int, itstep::Int, ntrac::Int)
    o.nsx = nsx
    o.nsy = nsy
    o.itstep = itstep
    o.ntrac = ntrac

    # optional parameters (default exists)
    o.subgrid = true
    o.subgrid_gaus = false
    o.score_method = "xcor"
    o.score_th0 = 0.8
    o.score_th1 = 0.7
    o.peak_inside_th = Float32(0.03) # unused if < 0
    o.min_contrast = Float32(-999.0) # unused if < 0
    o.vxch = -999.0 # unused if < 0
    o.vych = -999.0 # unused if < 0
    o.use_init_temp = false
end

"""
    set_optional!(o, subgrid, subgrid_gaus, score_method, score_th0, score_th1, peak_inside_th, min_contrast, vxch, vych, use_init_temp)

Sets optional tracking parameters.

# Arguments
- `o::VTT`: The object.
- `subgrid::Bool`: Whether to conduct subgrid tracking.
- `subgrid_gaus::Bool`: Whether subgrid peak finding is by gaussian.
- `score_method::String`: Scoring method (such as xcor for cross-correlation).
- `score_th0::AbstractFloat`: (Result screening parameter) Minimum score required for the first-time tracking.
- `score_th1::AbstractFloat`: (Result screening parameter) Minimum score required for the subsequent tracking.
- `peak_inside_th::Real`: An initial template is used only when
    it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.
- `min_contrast::Real`: An initial template is used only when 
    it has a difference in max and min greater than its value.
- `vxch::Float64`: (Result screening parameter) If positive, tracking result is rejected if the vx 
    chnages along trajecty greather than this value (thus used only when ntrac>=2). As a special case, 
    if the result of the second tracking is rejected, the first one is also rejected, since there is 
    no consecutive consistent result in this case.
- `vych::Float64`: (Result screening parameter) As vxch but for the y-component.
- `min_visible::Int`: Minimum number of visible values to calculate score when `chk_mask` is true.
"""
function set_optional!(o::VTT, subgrid::Bool, subgrid_gaus::Bool, score_method::String, score_th0::AbstractFloat, score_th1::AbstractFloat, peak_inside_th::Real, min_contrast::Real, vxch::Float64, vych::Float64, use_init_temp::Bool, min_visible::Int)
    o.subgrid = subgrid
    o.subgrid_gaus = subgrid_gaus
    o.score_method = score_method
    o.score_th0 = score_th0
    o.score_th1 = score_th1
    if peak_inside_th > 0.0
        o.chk_peak_inside = true
    else
        o.chk_peak_inside = false
    end
    o.peak_inside_th = Float32(peak_inside_th)

    if min_contrast > 0.0
        o.chk_min_contrast = true
    else
        o.chk_min_contrast = false
    end
    o.min_contrast = Float32(min_contrast)

    o.vxch = vxch # unused if < 0
    o.vych = vych # unused if < 0
    o.use_init_temp = use_init_temp
    o.min_visible = max(min_visible, 1)
end

"""To check whether a time index is valid. Returns `false` if valid, `true` if not."""
function inspect_t_index(o::VTT, tid::Int)
    stat = !(tid >= 1 && tid <= o.nt)
    return stat
end

"""
    get_zsub(o, tid, xi, yi)

Read out a template subimage from the image at `tid`.

The sub-image positons are specified at its center. (If the sub-image
size is even, with one more pix on the "left" / "bottom", by starting
from the index `xi-nsx/2`, `yi-nsy/2`).

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`, 
no data missing), `true` if not.
- `zsub::Matrix{Float32}`: Subimage at (x,y) = (xi, yi). 

# See Also
* [`get_zsub_view`](@ref)
"""
function get_zsub(o::VTT, tid::Int, xi::Int, yi::Int)
    stat = false
    nsx2, nsy2 = div(o.nsx,2), div(o.nsy,2)
    xi0, yi0 = xi - nsx2, yi - nsy2
    if xi0 < 1 || xi0 + o.nsx-1 > o.nx || yi0 < 1 || yi0 + o.nsy-1 > o.ny
        stat = true # sub-image is not within the original image
        return stat, nothing
    end
    zs = @inbounds o.z[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]
    if o.chk_zmiss
        stat = o.zmiss in zs
        if stat
            return stat, nothing
        end
    end
    return stat, zs
end

"""
    get_zsub_view(o, tid, xi, yi)

Like `get_zsub`, but returns a view.

# Notes
* This function returns a view of subimage. `get_zsub` returns a copy of subimage.

# See Also
* [`get_zsub`](@ref)
"""
function get_zsub_view(o::VTT, tid::Int, xi::Int, yi::Int)
    nsx2, nsy2 = div(o.nsx,2), div(o.nsy,2)
    xi0, yi0 = xi - nsx2, yi - nsy2
    if xi0 < 1 || xi0 + o.nsx-1 > o.nx || yi0 < 1 || yi0 + o.nsy-1 > o.ny
        return true, nothing # sub-image is not within the original image
    end
    zs = @inbounds @view o.z[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]
    if o.chk_zmiss
        if o.zmiss in zs
            return true, nothing
        end
    end
    return false, zs
end

"""
    get_zsub_visible(o, tid, xi, yi)

Read out a template subimage and subvisible from the image and visible at `tid`.

The sub-image positons are specified at its center. (If the sub-image
size is even, with one more pix on the "left" / "bottom", by starting
from the index `xi-nsx/2`, `yi-nsy/2`).

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`, 
no data missing), `true` if not.
- `zsub::Matrix{Float32}`: Subimage at (x,y) = (xi, yi). 

# See Also
* [`get_zsub_visible_view`](@ref)
"""
function get_zsub_visible(o::VTT, tid::Int, xi::Int, yi::Int)
    nsx2, nsy2 = div(o.nsx,2), div(o.nsy,2)
    xi0, yi0 = xi - nsx2, yi - nsy2
    if xi0 < 1 || xi0 + o.nsx-1 > o.nx || yi0 < 1 || yi0 + o.nsy-1 > o.ny
        return true, nothing, nothing # sub-image is not within the original image
    end
    zs = @inbounds o.z[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]
    if o.chk_zmiss
        if o.zmiss in zs
            return true, nothing, nothing
        end
    end

    visible = @inbounds o.visible[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]
    return false, zs, visible
end

"""
    get_zsub_visible_view(o, tid, xi, yi)

Like `get_zsub_visible`, but returns a view.

# Notes
* This function returns a view of zsub and visible. `get_zsub_visible` returns a copy of zsub and visible.

# See Also
* [`get_zsub_visible`](@ref)
"""
function get_zsub_visible_view(o::VTT, tid::Int, xi::Int, yi::Int)
    nsx2, nsy2 = div(o.nsx,2), div(o.nsy,2)
    xi0, yi0 = xi - nsx2, yi - nsy2
    if xi0 < 1 || xi0 + o.nsx-1 > o.nx || yi0 < 1 || yi0 + o.nsy-1 > o.ny
        return true, nothing, nothing # sub-image is not within the original image
    end

    zs = @inbounds @view o.z[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]
    if o.chk_zmiss
        if o.zmiss in zs
            return true, nothing, nothing
        end
    end

    visible = @inbounds @view o.visible[tid, yi0:yi0+o.nsy-1, xi0:xi0+o.nsx-1]

    return false, zs, visible
end

"""round to Int like C/C++"""
function roundInt(x::Real)
    return round(Int, x, RoundNearestTiesAway)
end

"""
    get_zsub_subgrid(o, tid, x, y)

Read out a template submimage from the image at `tid`.
Possibly at subgrid: Linearly interpolated, if
x or y has deviation from integer (bilinear if x and y).
Efficient: no unnecessary read-out is made.

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`, 
no data missing), `true` if not.
- `zsubg::Matrix{Float32}`: Subimage.
"""
function get_zsub_subgrid(o::VTT, tid::Int, x::Float64, y::Float64)
    xi, yi = roundInt(x), roundInt(y)
    dx, dy = x - xi, y - yi

    stat, zs = get_zsub_view(o, tid, xi, yi)
    if stat || (dx == 0.0 && dy == 0.0) # just on the grid
        return stat, zs
    end
    
    isx = Int(sign(dx))
    dx0 = Float32(abs(dx))
    dx1 = Float32(1.0)-dx0
    isy = Int(sign(dy))
    dy0 = Float32(abs(dy))
    dy1 = Float32(1.0)-dy0

    zsubg = zs * (dx1*dy1)
    if isx != 0
        stat, zsw1 = get_zsub_view(o, tid, xi+isx, yi)
        if stat
            return stat, nothing
        end
        zsubg .+= zsw1 * (dx0*dy1)
        if isy != 0
            stat, zsw1 = get_zsub_view(o, tid, xi+isx, yi+isy)
            if stat
                return stat, nothing
            end
            zsubg .+= zsw1 * (dx0*dy0)
        end
    end
    if isy != 0
        stat, zsw1 = get_zsub_view(o, tid, xi, yi+isy)
        if stat
            return stat, nothing
        end
        zsubg .+= zsw1 * (dx1*dy0)
    end
    return stat, zsubg
end
"""
    get_zsub_visible_subgrid(o, tid, x, y)

Read out a template submimage from the image at `tid`.
Possibly at subgrid: Linearly interpolated, if
x or y has deviation from integer (bilinear if x and y).
Efficient: no unnecessary read-out is made.

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`, 
no data missing), `true` if not.
- `zsubg::Matrix{Float32}`: Subimage.
"""
function get_zsub_visible_subgrid(o::VTT, tid::Int, x::Float64, y::Float64)
    xi, yi = roundInt(x), roundInt(y)
    dx, dy = x - xi, y - yi

    stat, zs, vissub = get_zsub_visible_view(o, tid, xi, yi)
    if stat || (dx == 0.0 && dy == 0.0) # just on the grid
        return stat, zs, vissub
    end

    isx = Int(sign(dx))
    dx0 = Float32(abs(dx))
    dx1 = Float32(1.0)-dx0
    isy = Int(sign(dy))
    dy0 = Float32(abs(dy))
    dy1 = Float32(1.0)-dy0

    zsubg = zs * (dx1*dy1)
    vissubg = vissub * (dx1*dy1)
    if isx != 0
        stat, zsw1, vissubw1 = get_zsub_visible_view(o, tid, xi+isx, yi)
        if stat
            return stat, nothing, nothing
        end
        zsubg .+= zsw1 * (dx0*dy1)
        vissubg .+= vissubw1 * (dx0*dy1)
        if isy != 0
            stat, zsw1, vissubw1 = get_zsub_visible_view(o, tid, xi+isx, yi+isy)
            if stat
                return stat, nothing, nothing
            end
            zsubg .+= zsw1 * (dx0*dy0)
            vissubg .+= vissubw1 * (dx0*dy0)
        end
    end
    if isy != 0
        stat, zsw1, vissubw1 = get_zsub_visible_view(o, tid, xi, yi+isy)
        if stat
            return stat, nothing, nothing
        end
        zsubg .+= zsw1 * (dx1*dy0)
        vissubg .+= vissubw1 * (dx1*dy0)
    end

    vissubg = Bool.(round.(Int8, vissubg))
    if !any(vissubg)
        return true, nothing, nothing
    end

    return stat, zsubg, vissubg
end

"""
    chk_zsub_peak_inside(o, zs)

Check whether the template subimage is peaked (maximized or minimized)
inside, and the peak is conspicuous enough, having a difference from the 
max or min on the sides greater than peak_inside_th*(inside_max - inside_min).

# Caution
If `o.peak_inside_th` < 0, no checking is conducted.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function chk_zsub_peak_inside(o::VTT, zs::AbstractMatrix{Float32})
    # find max and min along sides
    side_max = maximum([zs[begin,:]; zs[end,:]; zs[begin+1:end-1,begin]; zs[begin+1:end-1,end]])
    side_min = minimum([zs[begin,:]; zs[end,:]; zs[begin+1:end-1,begin]; zs[begin+1:end-1,end]])

    # find max and min inside the sides
    inner_max = maximum([zs[begin+1,begin+1:end-1]; zs[end-1,begin+1:end-1]; zs[begin+2:end-2,begin+1]; zs[begin+2:end-2,end-1]])
    inner_min = minimum([zs[begin+1,begin+1:end-1]; zs[end-1,begin+1:end-1]; zs[begin+2:end-2,begin+1]; zs[begin+2:end-2,end-1]])

    if (inner_max > side_max + o.peak_inside_th*(inner_max-inner_min) || inner_min < side_min - o.peak_inside_th*(inner_max-inner_min))
        return false # OK, because the max or min is inside and the difference from the max or min on the sides is not too tiny
    end
    return true
end

"""
    chk_zmiss_region(o, tid, k0, k1, l0, l1)

Check if there is data missing in the specified region at `tid`.

# Returns
- `stat::Bool`: `false` if there is no data missing, `true` if not.
"""
function chk_zmiss_region(o::VTT, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    return o.zmiss in @view o.z[tid, l0:l1, k0:k1]
end

"""
    sliding_xcor(o, sigx, xd, tid, k0, k1, l0, l1)

Sliding cross-correlation between the sugimage and image at `tid`.


# Returns
- `stat::Bool`: `false` if all the relevant data and regions are valid, so all the scores
    (xcor) are defined at all tested center locations; `true` if not.
- `scr::Matrix{Float32}`: Score array.
"""
function sliding_xcor(o::VTT, sigx::Real, xd::Matrix{Float32}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    nsx, nsy = o.nsx, o.nsy
    nsx2, nsy2 = div(nsx,2), div(nsy,2)
    nk = k1 - k0 + 1
    nl = l1 - l0 + 1
    nsxy = nsx * nsy
    k0 = k0 - nsx2
    l0 = l0 - nsy2
    scr = zeros(Float32, nl, nk)
    stat = ( k0 < 1 || k1+nsx2 > o.nx || l0 < 1 || l1+nsy2 > o.ny )
    if stat
        return stat, nothing
    end
    if o.chk_zmiss
        stat = chk_zmiss_region(o, tid, k0, k1+nsx2, l0, l1+nsy2)
        if stat
            return stat, nothing
        end
    end
    for l = 0:nl-1
        k = 0
        sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
        ymean = mean(sub_at_kl)
        yd = sub_at_kl .- ymean
        yysum = sum(yd.^2)
        xysum = sum(xd .* yd)
        vyy = yysum/nsxy
        vxy = xysum/nsxy
        scr[l+1,k+1] = vxy/sqrt(vyy)/sigx # cross-correlataion coef
        
        for k = 1:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            vyy = vyy + ymean^2 # mean(y^2) for previous k

            y_left = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k-1]
            y_right = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k+nsx-1]

            ydiff = sum(y_right) - sum(y_left)
            yydiff = sum(y_right.^2) - sum(y_left.^2)

            ymean += ydiff/nsxy # ymean is renewed.
            vyy = vyy + yydiff/nsxy - ymean^2 #new mean(y^2) - new ymean^2
            
            yd = sub_at_kl .- ymean
            xysum = sum(xd .* yd)
            vxy = xysum/nsxy
            scr[l+1,k+1] = vxy/sqrt(vyy)/sigx # cross-correlataion coef
        end
    end
    return stat, scr
end

"""
    sliding_ncov(o, sigx, xd, tid, k0, k1, l0, l1)

Sliding normalized covariance between the sugimage and image at `tid`.

Normalization is done by the sigma of the fist image : cov(x',y')/sigx^2
(in contrast to cov(x',y')/sigx/sigy in the correlation coefficient).

# Returns
- `stat::Bool`: `false` if all the relevant data and regions are valid, so all the scores
    (ncov) are defined at all tested center locations; `true` if not.
- `scr::Matrix{Float32}`: Score array.
"""
function sliding_ncov(o::VTT, sigx::Real, xd::Matrix{Float32}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    nsx, nsy = o.nsx, o.nsy
    nsx2, nsy2 = div(nsx,2), div(nsy,2)
    nk = k1 - k0 + 1
    nl = l1 - l0 + 1
    nsxy = nsx * nsy
    sigx2 = sigx^2
    k0 = k0 - nsx2
    l0 = l0 - nsy2
    scr = zeros(Float32, nl, nk)
    stat = ( k0 < 1 || k1+nsx2 > o.nx || l0 < 1 || l1+nsy2 > o.ny )
    if stat
        return stat, nothing
    end
    if o.chk_zmiss
        stat = chk_zmiss_region(o, tid, k0, k1+nsx2, l0, l1+nsy2)
        if stat
            return stat, nothing
        end
    end
    for l = 0:nl-1
        k = 0
        sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
        ymean = mean(sub_at_kl)
        yd = sub_at_kl .- ymean
        xysum = sum(xd .* yd)
        vxy = xysum/nsxy
        scr[l+1,k+1] = vxy/sigx2

        for k = 1:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]

            y_left = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k-1]
            y_right = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k+nsx-1]
            ydiff = sum(y_right) - sum(y_left)

            ymean = ymean + ydiff/nsxy # ymean is renewed.
            yd = sub_at_kl .- ymean
            xysum = sum(xd .* yd)
            vxy = xysum/nsxy
            scr[l+1,k+1] = vxy/sigx2
        end
    end
    return stat, scr
end

"""
    get_score_xcor(o, x, tid, k0, k1, l0, l1)

Conduct template matching, scoring by cross-correlation.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function get_score_xcor(o::VTT, x::AbstractMatrix{Float32}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    xm = mean(x)
    xd = x .- xm
    sigx = stdm(x, xm, corrected=false)
    stat, scr = sliding_xcor(o, sigx, xd, tid, k0, k1, l0, l1)    
    return stat, scr
end

"""
    get_score_ncov(o, x, tid, k0, k1, l0, l1)

Conduct template matching, scoring by normalized covariance.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function get_score_ncov(o::VTT, x::AbstractMatrix{Float32}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    xm = mean(x)
    xd = x .- xm
    sigx = stdm(x, xm, corrected=false)
    stat, scr = sliding_ncov(o, sigx, xd, tid, k0, k1, l0, l1)
    return stat, scr
end

"""
    get_score_xcor_with_visible(o, x, tid, k0, k1, l0, l1)

Conduct template matching, scoring by cross-correlation.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function get_score_xcor_with_visible(o::VTT, x::AbstractMatrix{Float32}, visible::Union{BitMatrix, AbstractMatrix{Bool}}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    nsx, nsy = o.nsx, o.nsy
    nsx2, nsy2 = div(nsx,2), div(nsy,2)
    nk = k1 - k0 + 1
    nl = l1 - l0 + 1
    k0 = k0 - nsx2
    l0 = l0 - nsy2
    scr = fill(Float32(o.fmiss), nl, nk)
    stat = ( k0 < 1 || k1+nsx2 > o.nx || l0 < 1 || l1+nsy2 > o.ny )
    if stat
        return stat, nothing
    end
    if o.chk_zmiss
        stat = chk_zmiss_region(o, tid, k0, k1+nsx2, l0, l1+nsy2)
        if stat
            return stat, nothing
        end
    end

    if all(@inbounds @view o.visible[tid, l0:l1+nsy2, k0:k1+nsx2])
        return get_score_xcor(o, x, tid, k0+nsx2, k1, l0+nsy2, l1)
    end

    allnan = true
    for l = 0:nl-1
        for k = 0:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_at_kl = @inbounds @view o.visible[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_and_visible = visible .* visible_at_kl
            if sum(visible_and_visible) < o.min_visible
                continue
            end
            scr[l+1,k+1] = cor(x[visible_and_visible], sub_at_kl[visible_and_visible])
            allnan = false
        end
    end

    if allnan
        return true, nothing
    end
    
    return stat, scr
end

"""
    get_score_ncov_with_visible(o, x, tid, k0, k1, l0, l1)

Conduct template matching, scoring by normalized covariance.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function get_score_ncov_with_visible(o::VTT, x::AbstractMatrix{Float32}, visible::Union{BitMatrix, AbstractMatrix{Bool}}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    nsx, nsy = o.nsx, o.nsy
    nsx2, nsy2 = div(nsx,2), div(nsy,2)
    nk = k1 - k0 + 1
    nl = l1 - l0 + 1
    k0 = k0 - nsx2
    l0 = l0 - nsy2
    scr = fill(Float32(o.fmiss), nl, nk)
    stat = ( k0 < 1 || k1+nsx2 > o.nx || l0 < 1 || l1+nsy2 > o.ny )
    if stat
        return stat, nothing
    end
    if o.chk_zmiss
        stat = chk_zmiss_region(o, tid, k0, k1+nsx2, l0, l1+nsy2)
        if stat
            return stat, nothing
        end
    end

    if all(@inbounds @view o.visible[tid, l0:l1+nsy2, k0:k1+nsx2])
        return get_score_ncov(o, x, tid, k0+nsx2, k1, l0+nsy2, l1)
    end

    allnan = true
    for l = 0:nl-1
        for k = 0:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_at_kl = @inbounds @view o.visible[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_and_visible = visible .* visible_at_kl
            if sum(visible_and_visible) < o.min_visible
                continue
            end
            x_valid = x[visible_and_visible]
            scr[l+1,k+1] = cov(x_valid, sub_at_kl[visible_and_visible], corrected=false)/std(x_valid, corrected=false)
            allnan = false
        end
    end

    if allnan
        return true, nothing
    end

    return stat, scr
end

"""Conduct template matching driver"""
function get_score(o::VTT, zs0::AbstractMatrix{Float32}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    if o.score_method == "xcor"
        stat, scr = get_score_xcor(o, zs0, tid, k0, k1, l0, l1)
    elseif o.score_method == "ncov"
        stat, scr = get_score_ncov(o, zs0, tid, k0, k1, l0, l1)
    else
        return true, nothing
    end
   return stat, scr
end

"""Conduct template matching driver with visible"""
function get_score_with_visible(o::VTT, zs0::AbstractMatrix{Float32}, visible::Union{BitMatrix, AbstractMatrix{Bool}}, tid::Int, k0::Int, k1::Int, l0::Int, l1::Int)
    if o.score_method == "xcor"
        stat, scr = get_score_xcor_with_visible(o, zs0, visible, tid, k0, k1, l0, l1)
    elseif o.score_method == "ncov"
        stat, scr = get_score_ncov_with_visible(o, zs0, visible, tid, k0, k1, l0, l1)
    else
        return true, nothing
    end
   return stat, scr
end

"""
    find_subgrid_peak_5pt_epara(c, l, r, b, t)

Find subgrid peak from 5 points with elliptic paraboloid.

Input c(enter), l(eft), r(ight), b(ottom), t(op) at
(0,0), (-1,0), (0,1), (0,-1), (0,1), respectively.
c must be greater than any of l,r,b,t.

Equation: z = -p(x-x0)^2 + -q(y-y0)^2 + r
            = -p x^2 + 2p x0 x - q y^2 + 2q y0 y + c

, where u = p x0, v = q y0, c = -p x0^2 - q y0^2 + r

# Returns
- `stat::Bool`: `false` if find peak successfully, `true` if not.
- `x0::float`: x of peak the location.
- `y0::float`: y of peak the location.

"""
function find_subgrid_peak_5pt_epara(c::Real, l::Real, r::Real, b::Real, t::Real)
    l = l-c
    r = r-c
    b = b-c
    t = t-c
    stat = !( l<=0.0 && r<=0.0 && b<=0.0 && t<=0.0 && ( l<0.0 || r<0.0 ) && ( b<0.0 || t<0.0) )
    if stat
        return stat, nothing, nothing, nothing
    end
    p = -(l+r)/2.0
    q = -(b+t)/2.0
    x0 = (r-l)/4p # --> |x0| < 0.5, if c >= [l,r,b,t] > 0
    y0 = (t-b)/4q # --> |y0| < 0.5, if c >= [l,r,b,t] > 0
    max_scr = c + p * x0^2 + q * y0^2 # r= c + p x0^2 + q y0^2
    return stat, x0, y0, max_scr
end

"""
    find_subgrid_peak_5pt_gaus(c, l, r, b, t)

Find subgrid peak from 5 points by interpolating with a 2D gaussian (for positive scores).

It is simply a log-version of the elliptic-paraboloid method
(`find_subgrid_peak_5pt_epara`). It appears that this method is
preferred in many PIVs over the elliptic-paraboloid method
(`find_subgrid_peak_5pt_epara`). In a test, there were not much
difference between them, though.

Input c(enter), l(eft), r(ight), b(ottom), t(op) at
(0,0), (-1,0), (0,1), (0,-1), (0,1), respectively.
all of them must be positive.
c must be greater than any of l,r,b,t.
"""
function find_subgrid_peak_5pt_gaus(c::Real, l::Real, r::Real, b::Real, t::Real)
    stat = !( c>0.0 && l>0.0 && r>0.0 && b>0.0 && t>0.0 ) # all must be >0
    if stat
        return stat, nothing, nothing, nothing
    end
    c = log(c)
    l = log(l) - c
    r = log(r) - c
    b = log(b) - c
    t = log(t) - c
    stat = !( l<=0.0 && r<=0.0 && b<=0.0 && t<=0.0 && ( l<0.0 || r<0.0 ) && ( b<0.0 || t<=0.0) )
    if stat
        return stat, nothing, nothing, nothing
    end
    p = -(l+r)/2.0
    q = -(b+t)/2.0
    x0 = (r-l)/4p # --> |x0| < 0.5, if c >= [l,r,b,t] > 0
    y0 = (t-b)/4q # --> |y0| < 0.5, if c >= [l,r,b,t] > 0
    max_scr = exp(c + p * x0^2 + q * y0^2)
    return stat, x0, y0, max_scr
end


"""print a 2D double array"""
function print_ary2d(a)
    for j = 1:size(a)[1]
        for i = 1:size(a)[2]
            @printf("%2.2f", a[j,i])
            print(" ")
        end
        print("\n")
    end
end


"""
    find_score_peak(o, scr, kw, lw)

Find the score peak and its location.

# Returns
- `stat::Bool`: `false` if the peak is inside; `true` if not.
- `kpi::Float64`: the peak location x.
- `lpi::Float64`: the peak location y.
- `scrp::Float64`: the peak score.
"""
function find_score_peak(o::VTT, scr::Matrix{Float32}, kw::Int, lw::Int)
    # find the max and its index
    if o.chk_mask
        l_and_k = findlast(x->x==maximum(filter(!isnan,scr)), scr)
    else
        l_and_k = findlast(x->x==maximum(scr), scr)
    end
    scrp = scr[l_and_k]
    lpi, kpi = l_and_k[1], l_and_k[2]

    # whether on the sides or not
    stat = ( kpi==1 || kpi==kw || lpi==1 || lpi==lw)
    if stat
        return stat, nothing, nothing, nothing
    end
    
    # subgrid determination
    if o.subgrid
        if o.subgrid_gaus
            stat, kp, lp, scrp = find_subgrid_peak_5pt_gaus(scr[lpi,kpi], scr[lpi,kpi-1], scr[lpi,kpi+1], scr[lpi-1,kpi], scr[lpi+1,kpi]) #[kl]p: relative to [kl]pi
        else
            stat, kp, lp, scrp = find_subgrid_peak_5pt_epara(scr[lpi,kpi], scr[lpi,kpi-1], scr[lpi,kpi+1], scr[lpi-1,kpi], scr[lpi+1,kpi]) #[kl]p: relative to [kl]pi
        end
        if stat
            return stat, nothing, nothing, nothing
        end
        kpi += kp
        lpi += lp
    end
    return stat, kpi, lpi, scrp
end


"""
    trac(o, tid, x, y[, vxg, vyg, out_subimage, out_score_ary])

Conduct tracking.

# Arguments
- `o::VTT`: The traking oject.
- `tid::Array{Integer,Any}`: Tracking initial time indices.
- `x::Array{Float64,Any}`: Tracking initial template-center x location (index-based; non-integer for subgrid).
- `y::Array{Float64,Any}`: Tracking initial template-center y location (index-based; non-integer for subgrid).
- `vxg::Array{Float64,Any}=nothing`: First guess of vx (to search around it). Can be 0.
- `vyg::Array{Float64,Any}=nothing`: First guess of vy (to search around it). Can be 0.
- `out_subimage::Bool=false`: Whether output subimages.
- `out_score_ary::Bool=false`: Whether output score arrays.
- `to_missing::Bool=true`: Whether output missing values as `missing`.

# Returns
- `count::Vector{Integer}`: [len] The number of successful tracking for each initial template.
- `tid::Matrix{Float64}`: [ntrac+1, len] time index of the trajectories (tid0 and subsequent ones).
- `x::Matrix{Float64}`: [ntrac+1, len] x locations of the trajectories (x0 and derived ones).
- `y::Matrix{Float64}`: [ntrac+1, len] y locations of trajectories (x0 and derived ones).
- `vx::Matrix{Float64}`: [ntrac, len] Derived x-velocity.
- `vy::Matrix{Float64}`: [ntrac, len] Derived y-velocity.
- `score::Matrix{Float64}`: [ntrac, len] Scores along the trajectory (max values, possibly at subgrid).
- `zss::Array{Float32,4}`: [nsx, nsy, ntrac+1, len] (optional, if non-`nothing`)
    (Diagnosis output if wanted) The subimages along the track.
- `score_arry::Array{Float64,4}`: [(x-sliding size, y-sliding size, ntrac+1, len] (optional, if non-`nothing`)
    (Diagnosis output if wanted) The entire scores.
"""
function trac(o::VTT, tid, x, y; vxg=nothing, vyg=nothing, out_subimage::Bool=false, out_score_ary::Bool=false, to_missing::Bool=true)
    !o.setuped && throw(ArgumentError("Need to call #setup in advance"))
    sh = size(x)
    if typeof(tid) === Int64
        tid = fill(tid, sh...)
    else
        size(tid) !== sh && throw(ArgumentError("Shape miss-match (x)"))
    end

    size(y) != sh && throw(ArgumentError("Shape miss-match (y)"))
    
    if isnothing(vxg)
        vxg = zeros(sh...)
    else
        size(vxg) !== sh && throw(ArgumentError("Shape miss-match (vxg)"))
    end

    if isnothing(vyg)
        vyg = zeros(sh...)
    else
        size(vyg) !== sh && throw(ArgumentError("Shape miss-match (vyg)"))
    end

    count, status, tid, x, y, vx, vy, score, zss, score_ary = do_tracking(o, vec(tid), vec(x), vec(y), vec(vxg), vec(vyg), out_subimage, out_score_ary)

    if to_missing
        fmiss = o.fmiss
        imiss = o.imiss
        tid = Array{Union{Missing, Int},2}(tid)
        tid[tid.==imiss] .= missing
        x = Array{Union{Missing, Float64},2}(x)
        x[x.==fmiss] .= missing
        y = Array{Union{Missing, Float64},2}(y)
        y[y.==fmiss] .= missing
        vx = Array{Union{Missing, Float64},2}(vx)
        vx[vx.==fmiss] .= missing
        vy = Array{Union{Missing, Float64},2}(vy)
        vy[vy.==fmiss] .= missing
        score = Array{Union{Missing, Float32},2}(score)
        score[score.==fmiss] .= missing
        zss = Array{Union{Missing, Float32},4}(zss)
        zss[zss.==o.zmiss] .= missing
        score_ary = Array{Union{Missing, Float32},4}(score_ary)
        score_ary[score_ary.==Float32(fmiss)] .= missing
    end

    if length(sh) >= 2
        # reshape outputs based on the shape of inputs
        count = reshape(count, sh...)
        status = reshape(status, sh...)
        tid = reshape(tid, size(tid)[1], sh...)
        x = reshape(x, size(x)[1], sh...)
        y = reshape(y, size(y)[1], sh...)
        vx = reshape(vx, size(vx)[1], sh...)
        vy = reshape(vy, size(vy)[1], sh...)
        score = reshape(score, size(score)[1], sh...)
        if out_subimage
            zss = reshape(zss, size(zss)[1:end-1]..., sh...)
        end
        if out_score_ary
            score_ary = reshape(score_ary, size(score_ary)[1:end-1]..., sh...)
        end
    end
    return count, status, tid, x, y, vx, vy, score, zss, score_ary
end


"""
    do_tracking(o, tid0, x0, y0, vx0, vy0, out_subimage, out_score_ary)

Conduct tracking (core).

# Arguments
- `o::VTT`: The traking oject.
- `tid0::Vector{Integer}`: Tracking initial time indices.
- `x0::Vector{Float64}`: Tracking initial template-center x location (index-based; non-integer for subgrid).
- `y0::Vector{Float64}`: Tracking initial template-center y location (index-based; non-integer for subgrid).
- `vx0g::Vector{Float64}`: First guess of vx (to search around it). Can be 0.
- `vy0g::Vector{Float64}`: First guess of vy (to search around it). Can be 0.
- `out_subimage::Bool`: Whether output subimages.
- `out_score_ary::Bool`: Whether output score arrays.

# Returns
- `count::Matrix{Float64}`: (len: len) The number of successful tracking for each initial template.
- `tid::Matrix{Float64}`: (len: (ntrac+1)*len) Time index of the trajectories (tid0 and subsequent ones).
- `x::Matrix{Float64}`: (len: (ntrac+1)*len) x locations of the trajectories (x0 and derived ones).
- `y::Matrix{Float64}`: (len: (ntrac+1)*len) y locations of trajectories (x0 and derived ones).
- `vx::Matrix{Float64}`: (len: ntrac*len) Derived x-velocity.
- `vy::Matrix{Float64}`: (len: ntrac*len) Derived y-velocity.
- `score::Matrix{Float64}`: (len: ntrac*len)  Scores along the trajectory (max values, possibly at subgrid).
- `zss::Array{Float32,4}`: (optional, if non-`nothing`) (Diagnosis output if wanted) The subimages along the track (1D pointer for 4D array; nsx * nsy * (ntrac+1) * len.
- `score_ary::Array{Float32,4}`: (optional, if non-`nothing`) (Diagnosis output if wanted) the entire scores (1D pointer for 4D array; (x-sliding size) * (y-sliding size) * (ntrac+1) * len.
"""
function do_tracking(o::VTT, tid0, x0, y0, vx0, vy0, out_subimage::Bool, out_score_ary::Bool)
    len = length(tid0)
    fmiss = o.fmiss
    imiss = o.imiss

    length(x0) != len && throw(ArgumentError("invalid x0 length"))
    length(y0) != len && throw(ArgumentError("invalid y0 length"))
    length(vx0) != len && throw(ArgumentError("invalid vx0 length"))
    length(vy0) != len && throw(ArgumentError("invalid vy0 length"))

    shape0 = [o.ntrac, len]
    shape1 = [o.ntrac+1, len]
    count = zeros(Int, len)
    tid = fill(imiss, shape1...)
    x = fill(fmiss, shape1...)
    y = fill(fmiss, shape1...)
    vx = fill(fmiss, shape0...)
    vy = fill(fmiss, shape0...)
    score = fill(Float32(fmiss), shape0...)
    zs0 = fill(o.zmiss, o.nsy, o.nsx)

    if out_subimage
        shape2 = [o.ntrac+1, o.nsy, o.nsx, len]
        zss = fill(o.zmiss, shape2...)
    else
        zss = nothing
    end

    if out_score_ary
        shape2 = [o.ntrac, 2o.iyhw+1, 2o.ixhw+1, len]
        score_ary = fill(Float32(fmiss), shape2...)
    else
        score_ary = nothing
    end

    ixhw = o.ixhw
    iyhw = o.iyhw
    kw = 2ixhw + 1
    lw = 2iyhw + 1
    t = o.t
    itstep = o.itstep
    chk_vchange = (o.vxch > 0.0 && o.vych > 0.0)

    status = zeros(len)

    # record initial data
    if o.subgrid # initial position
        x[1,:] .= x0[:]
        y[1,:] .= y0[:]
    else
        x[1,:] .= roundInt.(x0[:])
        y[1,:] .= roundInt.(y0[:])
    end
    tid[1,:] .= tid0[:]

    for m = 1:len
        for j = 1:o.ntrac
            if status[m] != 0
                continue
            end
        
            xcur = @inbounds x[j,m] # x current
            ycur = @inbounds y[j,m] # y current

            # record initial data
            tidf = tid0[m] + (j-1)*itstep  # index of the tracking start time
            tidl = tidf + itstep      # index of the tracking end time
            stat = inspect_t_index(o, tidf)
            if stat
                status[m] = 1
                # @info "(m=$m) Stop tracking at checkpoint 1 (during `inspect_t_index` of `tidf`)"
                continue
            end
            if j == 1 || !o.use_init_temp 
                if o.chk_mask
                    stat, zs0, visible = get_zsub_visible_subgrid(o, tidf, xcur, ycur)
                    if !stat
                        stat = sum(visible) < o.min_visible
                    end
                else
                    stat, zs0 = get_zsub_subgrid(o, tidf, xcur, ycur)
                end
            end

            if stat
                status[m] = 2
                # @info "(m=$m) Stop tracking at checkpoint 2 (during `get_zsub_subgrid`)"
                continue
            end

            if o.chk_min_contrast
                if maximum(zs0) - minimum(zs0) < o.min_contrast
                    status[m] = 3
                    # @info "(m=$m) Stop tracking at checkpoint 3 (during `min_contrast` check)"
                    continue
                end
            end

            if o.chk_peak_inside
                stat = chk_zsub_peak_inside(o, zs0)
                if stat
                    status[m] = 4
                    # @info "(m=$m) Stop tracking at checkpoint 4 (during `chk_zsub_peak_inside`)"
                    continue
                end
            end

            if out_subimage
                if j == 1 || !o.use_init_temp
                    zss[j,:,:,m] = zs0
                else
                    stat, zsw = get_zsub_subgrid(o, tidf, xcur, ycur)
                    zss[j,:,:,m] = zsw
                end
            end

            # inspect the tracking end time
            stat = inspect_t_index(o, tidl)
            if stat
                status[m] = 5
                # @info "(m=$m) Stop tracking at checkpoint 5 (during `inspect_t_index` of `tidl`)"
                continue
            end
            dt = t[tidl] - t[tidf] # time diff. can be negative
            if j == 1
                vxg = vx0[m]
                vyg = vy0[m]
            else 
                vxg = vx[j-1,m] # previous result
                vyg = vy[j-1,m] # previous result
            end
            kc = roundInt(xcur + vxg*dt)
            lc = roundInt(ycur + vyg*dt)
            
            if o.chk_mask
                stat, scr = get_score_with_visible(o, zs0, visible, tidl, kc-ixhw, kc+ixhw, lc-iyhw, lc+iyhw)
            else
                stat, scr = get_score(o, zs0, tidl, kc-ixhw, kc+ixhw, lc-iyhw, lc+iyhw)
            end

            if stat
                status[m] = 6
                # @info "(m=$m) Stop tracking at checkpoint 6 (during `get_score`)"
                continue
            end

            if out_score_ary
                score_ary[j,:,:,m] .= scr
            end

            # print_dary2d("**score", "%7.3f", scr, kw, lw );
            stat, xp, yp, sp = find_score_peak(o, scr, kw, lw)
            if stat
                status[m] = 7
                # @info "(m=$m) Stop tracking at checkpoint 7 (during `find_score_peak`)"
                continue
            end
            if ((j==1 && sp<o.score_th0) || (j>1 && sp<o.score_th1))
                status[m] = 8
                # @info "(m=$m) Stop tracking at checkpoint 8 (during `score_th0` or `score_th1` check)"
                continue
            end
            score[j,m] = Float32(sp)
            xw = xp + kc - 1 - ixhw # next position (x-axis)
            yw = yp + lc - 1 - iyhw # next position (y-axis)
            vxw = (xw - x[j,m])/dt # velocity (x-axis)
            vyw = (yw - y[j,m])/dt # velocity (y-axis)
            if chk_vchange && j > 1
                stat = (abs(vxw - vx[j-1,m]) > o.vxch || abs(vyw - vy[j-1,m]) > o.vych)
                if stat
                    if j == 2 # invalidate the j==1 results too
                        count[m] = 0
                        x[1,m] = y[1,m] = o.fmiss
                        vx[1,m] = vy[1,m] = o.fmiss
                    end
                    status[m] = 9
                    # @info "(m=$m) Stop tracking at checkpoint 9 (during `chk_vchange`)"
                    continue
                end
            end
            count[m] = j + 1 # of valid tracking
            tid[j+1,m] = tidl
            x[j+1,m] = xw
            y[j+1,m] = yw
            vx[j,m] = vxw
            vy[j,m] = vyw
            if out_subimage && j == o.ntrac # last sub image
                xcur = @inbounds x[j+1,m]
                ycur = @inbounds y[j+1,m]
                stat, zs0 = get_zsub_subgrid(o, tidf, xcur, ycur)
                zss[j+1,:,:,m] = zs0
            end
        end
    end   
    return count, status, tid, x, y, vx, vy, score, zss, score_ary
end
end
