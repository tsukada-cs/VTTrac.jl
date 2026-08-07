module VTTrac
export VTT, setup, trac, set_ixyhw_from_v, set_ixyhw_directly

using Statistics


"""
    VTT(z; t=nothing, mask=nothing, zmiss=nothing, fmiss=-999.0, imiss=-999)

Sets data for tracking (you need to set parameters separately with [`setup`](@ref)).

# Arguments
- `z::AbstractArray{Float32,3}`: Array of image-like data (in dimensions [time, y, x]), with at least 2 time steps. `z[i]` contains `i`-th image data.
- `t::Union{AbstractVector{Float64}, Nothing}=nothing`: Times at which the images are for. Defaults to `0:size(z,1)-1`.
- `mask::Union{Array{Bool,3}, Nothing}=nothing`: Mask to ignore when calculating score (true positions are ignored).
- `zmiss::Union{Real, Nothing}=nothing`: Missing value used in `z`.
- `fmiss::Real=-999.0`: Missing value to be set for Real.
- `imiss::Integer=-999`: Missing value to be set for Integer.
"""
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
    Sth0::Float64
    Sth1::Float64
    chk_peak_inside::Bool
    peak_inside_th::Float32 #threshold for the peak-inside screening(unused if<0)
    chk_Cth::Bool
    Cth::Float32 # minimum contrast in the template (unused if=<0)
    use_init_temp::Bool # if true, always use initial template submimages
    setuped::Bool
    min_samples::Int # minimum number of visible values to calculate score when `chk_mask` is true.

    function VTT(z::AbstractArray{Float32,3}; t::Union{AbstractVector{Float64}, Nothing}=nothing, mask::Union{BitArray{3}, Array{Bool,3}, Nothing}=nothing, zmiss::Union{Real, Nothing}=nothing, fmiss::Real=-999.0, imiss::Int=-999)
        o = new()
        o.z = z
        o.nt, o.ny, o.nx = size(z)
        o.nt < 2 && throw(ArgumentError("`z` must have at least 2 time steps (size(z, 1) >= 2)"))
        if isnothing(t)
            t = Vector{Float64}([0:o.nt-1;])
        end
        size(t) != (o.nt,) && throw(ArgumentError("`size(t)` must be `(size(z)[begin],1)`"))
        o.t = t

        if isnothing(zmiss)
            o.chk_zmiss = false
            o.zmiss = Float32(0.0)
        else
            o.chk_zmiss = true
            o.zmiss = Float32(zmiss)
        end
        o.fmiss = fmiss
        o.imiss = imiss

        if isnothing(mask)
            o.chk_mask = false
        else
            size(mask) != size(z) && throw(ArgumentError("`size(mask)` must be `size(z)`"))
            if any(mask)
                o.chk_mask = true
                o.visible = .!mask
            else
                o.chk_mask = false
            end
        end

        o.setuped = false
        return o
    end
end

"""
    setup(o, nsx, nsy, vxhw, vyhw, [ixhw, iyhw, subgrid, subgrid_gaus,
        itstep, ntrac, score_method, Sth0, Sth1, vxch, vych,
        peak_inside, peak_inside_th, Cth, use_init_temp, min_samples])

Setup for tracking.

# Arguments
- `o::VTT`: The object.
- `nsx::Integer, nsy::Integer`: Subimage x & y sizes (x:1st, y:2nd dim).
- `vxhw::Union{Real, Nothing}, vyhw::Union{Real, Nothing}`: (either `v[xy]hw` or `i[xy]hw` are MANDATORY).
    the dimensions along which to perform the computation.
    search velocity range half sizes to set `i[xy]hw`.
    Seach at least to cover +-v?hw around the first guess or previous step.
    (the result can be outside the range.)
- `ixhw::Union{Integer, Nothing}, iyhw::Union{Integer, Nothing}`: (either `v[xy]hw` or `i[xy]hw` are MANDATORY)
    Max displacement fro template match (can be set indirecly through `v[xy]hw`).
- `subgrid::Bool=true`: Whether to conduct subgrid tracking.
- `subgrid_gaus::Bool=true`: Whether subgrid peak finding is by gaussian.
- `itstep::Integer=1`: Step of `t`'s used (skip if >1).
- `ntrac::Integer=2`: Max tracking times from initial loc.
- `score_method::String="xcor"`: `"xcor"` for cross-correlation, `"ncov"` for normalized covariance.
- `Sth0::Real=0.8`: The minimum score required for the 1st tracking.
- `Sth1::Real=0.7`: The minimum score required for subsequent tracking.
- `vxch::Union{Real, Nothing}=nothing`: If non-`nothing`, the max tolerant vx
    change between two consecutive tracking. Screening by velocity change is only
    active when both `vxch` and `vych` are set (both non-`nothing` and positive).
- `vych::Union{Real, Nothing}=nothing`: If non-`nothing`, the max tolerant vy
    change between two consecutive tracking. See `vxch`.
- `peak_inside_th::Union{Real, Nothing}=nothing`: If non-`nothing`, an initial template is used only when
    it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.
- `Cth::Union{Real, Nothing}=nothing`: If non-`nothing`, an initial template is used only when 
    it has a difference in max and min greater than its value.
- `min_samples::Integer=1`: Minimum number of visible values to calculate score when `chk_mask` is true.
"""
function setup(
    o::VTT, nsx::Integer, nsy::Integer;
    vxhw::Union{Real, Nothing}=nothing, vyhw::Union{Real, Nothing}=nothing,
    ixhw::Union{Integer, Nothing}=nothing, iyhw::Union{Integer, Nothing}=nothing,
    subgrid::Bool=true, subgrid_gaus::Bool=false,
    itstep::Integer=1, ntrac::Integer=2,
    score_method::String="xcor",
    Sth0::Real=0.8, Sth1::Real=0.7,
    vxch::Union{Real, Nothing}=nothing, vych::Union{Real, Nothing}=nothing,
    peak_inside_th::Union{Real, Nothing}=nothing,
    Cth::Union{Real, Nothing}=nothing,
    use_init_temp::Bool=false,
    min_samples::Integer=1
    )
    nsx, nsy = Int(nsx), Int(nsy)
    itstep, ntrac, min_samples = Int(itstep), Int(ntrac), Int(min_samples)
    Sth0, Sth1 = Float64(Sth0), Float64(Sth1)

    (nsx < 1 || nsy < 1) && throw(ArgumentError("`nsx` and `nsy` must be positive"))
    score_method ∉ ("xcor", "ncov") && throw(ArgumentError("`score_method` must be either \"xcor\" or \"ncov\""))

    o.nsx = nsx
    o.nsy = nsy
    o.dtmean = itstep * (o.t[end]-o.t[begin]) / (o.nt-1)
    if vxhw !== nothing
        ixhw !== nothing && throw(ArgumentError("`v[xy]hw` and `i[xy]hw` must not be set simultaneously"))
        vyhw === nothing && throw(ArgumentError("vxhw and vyhw must be set simultaneously"))
        set_ixyhw_from_v(o, vxhw, vyhw)
    elseif ixhw !== nothing
        iyhw === nothing && throw(ArgumentError("ixhw and iyhw must be set simultaneously"))
        set_ixyhw_directly(o, ixhw, iyhw)
    else
        throw(ArgumentError("either `i[xy]hw` or `v[xy]hw` must be specified"))
    end

    vxch === nothing && (vxch = -999.0) # <=0 for nothing (not to set)
    vych === nothing && (vych = -999.0) # <=0 for nothing (not to set)

    if isnothing(peak_inside_th)
        peak_inside_th = -1.0  # negative, meaning unused
    end
    peak_inside_th = Float32(peak_inside_th)
    if peak_inside_th > 0.0 && (nsx < 3 || nsy < 3)
        throw(ArgumentError("`nsx` and `nsy` must be at least 3 when `peak_inside_th` is set, since the peak-inside check needs an interior region"))
    end
    if isnothing(Cth)
        Cth = -1.0   # negative, meaning unused
    end
    Cth = Float32(Cth)

    set_basic!(o, nsx, nsy, itstep, ntrac)
    set_optional!(o, subgrid, subgrid_gaus, score_method, Sth0, Sth1, peak_inside_th, Cth, vxch, vych, use_init_temp, min_samples)
    o.setuped = true
end

"""
    set_ixyhw_from_v(o, vxhw, vyhw)

Sets the tracking parameters `i[xy]hw` from velocities (v[xy]hw).

# Arguments
- `o::VTT`: The object.
- `vxhw::Real`: The range over which vx is searched around initial guess.
- `vyhw::Real`: The range over which vy is searched around initial guess.
"""
function set_ixyhw_from_v(o::VTT, vxhw::Real, vyhw::Real)
    o.vxhw = Float64(vxhw)
    o.vyhw = Float64(vyhw)
    o.ixhw = ceil(Int, abs(o.vxhw * o.dtmean)) + 1 # max displacement
    o.iyhw = ceil(Int, abs(o.vyhw * o.dtmean)) + 1 # +1 is margin to find peak
end

"""
    set_ixyhw_directly(o, ixhw, iyhw)

Sets the tracking parameters `i[xy]hw`.

# Arguments
- `o::VTT`: The object.
- `ixhw::Integer`: The range over which next x is searched around initial guess.
- `iyhw::Integer`: The range over which next y is searched around initial guess.
"""
function set_ixyhw_directly(o::VTT, ixhw::Integer, iyhw::Integer)
    o.ixhw = Int(ixhw)
    o.iyhw = Int(iyhw)
    o.vxhw = (o.ixhw - 1) / o.dtmean # max displacement
    o.vyhw = (o.iyhw - 1) / o.dtmean # -1 is from margin to find peak
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
    o.Sth0 = 0.8
    o.Sth1 = 0.7
    o.peak_inside_th = Float32(0.03) # unused if < 0
    o.Cth = Float32(-999.0) # unused if < 0
    o.vxch = -999.0 # unused if < 0
    o.vych = -999.0 # unused if < 0
    o.use_init_temp = false
end

"""
    set_optional!(o, subgrid, subgrid_gaus, score_method, Sth0, Sth1, peak_inside_th, Cth, vxch, vych, use_init_temp)

Sets optional tracking parameters.

# Arguments
- `o::VTT`: The object.
- `subgrid::Bool`: Whether to conduct subgrid tracking.
- `subgrid_gaus::Bool`: Whether subgrid peak finding is by gaussian.
- `score_method::String`: Scoring method (such as xcor for cross-correlation).
- `Sth0::Float64`: (Result screening parameter) Minimum score required for the first-time tracking.
- `Sth1::Float64`: (Result screening parameter) Minimum score required for the subsequent tracking.
- `peak_inside_th::Real`: An initial template is used only when
    it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.
- `Cth::Real`: An initial template is used only when 
    it has a difference in max and min greater than its value.
- `vxch::Float64`: (Result screening parameter) If positive, tracking result is rejected if the vx
    changes along trajectory greater than this value (thus used only when ntrac>=2). As a special case,
    if the result of the second tracking is rejected, the first one is also rejected, since there is
    no consecutive consistent result in this case.
- `vych::Float64`: (Result screening parameter) As vxch but for the y-component.
- `min_samples::Int`: Minimum number of visible values to calculate score when `chk_mask` is true.
"""
function set_optional!(o::VTT, subgrid::Bool, subgrid_gaus::Bool, score_method::String, Sth0::Float64, Sth1::Float64, peak_inside_th::Real, Cth::Real, vxch::Float64, vych::Float64, use_init_temp::Bool, min_samples::Int)
    o.subgrid = subgrid
    o.subgrid_gaus = subgrid_gaus
    o.score_method = score_method
    o.Sth0 = Sth0
    o.Sth1 = Sth1
    if peak_inside_th > 0.0
        o.chk_peak_inside = true
    else
        o.chk_peak_inside = false
    end
    o.peak_inside_th = Float32(peak_inside_th)

    if Cth > 0.0
        o.chk_Cth = true
    else
        o.chk_Cth = false
    end
    o.Cth = Float32(Cth)

    o.vxch = vxch # unused if < 0
    o.vych = vych # unused if < 0
    o.use_init_temp = use_init_temp
    o.min_samples = max(min_samples, 1)
end

"""To check whether a time index is valid. Returns `false` if valid, `true` if not."""
function inspect_t_index(o::VTT, tid::Int)
    stat = !(tid >= 1 && tid <= o.nt)
    return stat
end

"""
    get_zsub_view(o, tid, xi, yi)

Read out a template subimage (as a view) from the image at `tid`.

The sub-image positons are specified at its center. (If the sub-image
size is even, with one more pix on the "left" / "bottom", by starting
from the index `xi-nsx/2`, `yi-nsy/2`).

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`,
no data missing), `true` if not.
- `zsub::Matrix{Float32}`: Subimage at (x,y) = (xi, yi).
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
    get_zsub_visible_view(o, tid, xi, yi)

Read out a template subimage and its visible mask (as views) from the image and
visible mask at `tid`.

The sub-image positons are specified at its center. (If the sub-image
size is even, with one more pix on the "left" / "bottom", by starting
from the index `xi-nsx/2`, `yi-nsy/2`).

# Returns
- `stats::Bool`: `false` if successful (specified region is valid and, if `chk_zmiss`,
no data missing), `true` if not.
- `zsub::Matrix{Float32}`: Subimage at (x,y) = (xi, yi).
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

Check if the subimage has a prominent peak or trough in its interior.
A peak or trough is considered prominent if its value differs from the boundary 
extremes by more than `o.peak_inside_th` times the interior data range.

# Caution
If `o.peak_inside_th` < 0, no checking is conducted.

# Returns
- `stat::Bool`: `false` if passed the check, `true` if not.
"""
function chk_zsub_peak_inside(o::VTT, zs::AbstractMatrix{Float32})
    # find max and min along sides
    side_max = max(maximum(@view zs[begin, :]), maximum(@view zs[end, :]),
                   maximum(@view zs[begin+1:end-1, begin]), maximum(@view zs[begin+1:end-1, end]))
    side_min = min(minimum(@view zs[begin, :]), minimum(@view zs[end, :]),
                   minimum(@view zs[begin+1:end-1, begin]), minimum(@view zs[begin+1:end-1, end]))

    # find max and min inside the sides
    inner_max = maximum(@view zs[begin+1:end-1,begin+1:end-1])
    inner_min = minimum(@view zs[begin+1:end-1,begin+1:end-1])

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
Sum of squared deviations of `sub` from `ymean`, and of the cross deviations between
`sub` and the (already demeaned) template `xd`, computed in a single allocation-free pass.

Demeaning `sub` before squaring/multiplying (rather than accumulating raw sums of squares)
avoids catastrophic cancellation for data with a large mean offset relative to its variance.
"""
@inline function xcor_sums(sub::AbstractMatrix{Float32}, xd::Matrix{Float32}, ymean::Float32)
    vyy_sum = 0.0f0
    xysum = 0.0f0
    @inbounds for i in eachindex(sub, xd)
        d = sub[i] - ymean
        vyy_sum += d*d
        xysum += xd[i]*d
    end
    return vyy_sum, xysum
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
    scr = zeros(Float32, nl, nk)
    for l = 0:nl-1
        k = 0
        sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
        ymean = mean(sub_at_kl)
        vyy_sum, xysum = xcor_sums(sub_at_kl, xd, ymean)
        vyy = vyy_sum/nsxy
        if vyy <= 0
            scr[l+1,k+1] = o.fmiss
        else
            vxy = xysum/nsxy
            scr[l+1,k+1] = vxy/sqrt(vyy)/sigx # cross-correlataion coef
        end

        for k = 1:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]

            y_left = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k-1]
            y_right = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k+nsx-1]

            ydiff = sum(y_right) - sum(y_left)
            ymean += ydiff/nsxy # ymean is renewed (running sum of differences; numerically stable).

            # `vyy` (variance of the current window) is recomputed directly from the
            # demeaned window each time, rather than incrementally from running sums of
            # squares. The latter (E[y^2] - E[y]^2) suffers catastrophic cancellation for
            # data with a large mean offset relative to its variance, which could previously
            # make `vyy` spuriously go negative and, once reset to 0, permanently corrupt
            # the running state for the remainder of the row.
            vyy_sum, xysum = xcor_sums(sub_at_kl, xd, ymean)
            vyy = vyy_sum/nsxy

            if vyy <= 0
                scr[l+1,k+1] = o.fmiss
                continue
            end

            vxy = xysum/nsxy
            scr[l+1,k+1] = vxy/sqrt(vyy)/sigx # cross-correlataion coef
        end
    end
    return stat, scr
end

"""
    sliding_ncov(o, sigx, xd, tid, k0, k1, l0, l1)

Sliding normalized covariance between the subimage and image at `tid`.

Normalization is done by the sigma of the first image : cov(x',y')/sigx^2
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
    scr = zeros(Float32, nl, nk)
    for l = 0:nl-1
        k = 0
        sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
        ymean = mean(sub_at_kl)
        _, xysum = xcor_sums(sub_at_kl, xd, ymean)
        vxy = xysum/nsxy
        scr[l+1,k+1] = vxy/sigx2

        for k = 1:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]

            y_left = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k-1]
            y_right = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k+nsx-1]
            ydiff = sum(y_right) - sum(y_left)

            ymean = ymean + ydiff/nsxy # ymean is renewed.
            _, xysum = xcor_sums(sub_at_kl, xd, ymean)
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
Sample size and the sums needed for masked `xcor`/`ncov` scores, restricted to positions
where both `visible` and `visible_at_kl` are true, computed in a two-pass (mean, then
demeaned sums), allocation-free way. Demeaning before squaring/multiplying (rather than
accumulating raw sums of squares) avoids catastrophic cancellation for data with a large
mean offset relative to its variance (see `xcor_sums`).

# Returns
- `n::Int`: Number of jointly-visible samples.
- `sxx, syy, sxy::Float32`: Sums of squared/cross deviations from the (masked) means.
"""
@inline function masked_sums(x::AbstractMatrix{Float32}, sub::AbstractMatrix{Float32},
                              visible::AbstractMatrix{Bool}, visible_at_kl::AbstractMatrix{Bool})
    n = 0
    sx = 0.0f0
    sy = 0.0f0
    @inbounds for i in eachindex(x, sub, visible, visible_at_kl)
        if visible[i] && visible_at_kl[i]
            n += 1
            sx += x[i]
            sy += sub[i]
        end
    end
    n == 0 && return 0, 0.0f0, 0.0f0, 0.0f0
    mx = sx/n
    my = sy/n
    sxx = 0.0f0
    syy = 0.0f0
    sxy = 0.0f0
    @inbounds for i in eachindex(x, sub, visible, visible_at_kl)
        if visible[i] && visible_at_kl[i]
            dx = x[i] - mx
            dy = sub[i] - my
            sxx += dx*dx
            syy += dy*dy
            sxy += dx*dy
        end
    end
    return n, sxx, syy, sxy
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

    scr = fill(Float32(o.fmiss), nl, nk)
    allnan = true
    for l = 0:nl-1
        for k = 0:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_at_kl = @inbounds @view o.visible[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            n, sxx, syy, sxy = masked_sums(x, sub_at_kl, visible, visible_at_kl)
            if n < o.min_samples
                continue
            end
            scr[l+1,k+1] = sxy / sqrt(sxx*syy) # cross-correlation coef; matches `cor(x,y)`
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

    scr = fill(Float32(o.fmiss), nl, nk)
    allnan = true
    for l = 0:nl-1
        for k = 0:nk-1
            sub_at_kl = @inbounds @view o.z[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            visible_at_kl = @inbounds @view o.visible[tid, l0+l:l0+l+nsy-1, k0+k:k0+k+nsx-1]
            n, sxx, _, sxy = masked_sums(x, sub_at_kl, visible, visible_at_kl)
            if n < o.min_samples
                continue
            end
            # cov(x,y)/var(x), consistent with `sliding_ncov` (matches
            # `cov(x,y,corrected=false)/var(x,corrected=false)`)
            scr[l+1,k+1] = sxy / sxx
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
(0,0), (-1,0), (1,0), (0,-1), (0,1), respectively.
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
(0,0), (-1,0), (1,0), (0,-1), (0,1), respectively.
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
    stat = !( l<=0.0 && r<=0.0 && b<=0.0 && t<=0.0 && ( l<0.0 || r<0.0 ) && ( b<0.0 || t<0.0) )
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


"""
    find_score_peak(o, scr, kw, lw)

Find the maximum score and its subgrid location in the score matrix `scr`.

If the maximum score is located at the boundary of the matrix, subgrid interpolation
cannot be performed. In this case, the function returns the integer coordinates of
the maximum and sets `stat` to `true`. For interior peaks, a 5-point Gaussian
interpolation is applied. If the interpolation is successful, `stat` is `false` and
the subgrid coordinates are returned. If interpolation fails, `stat` is set to `true`.


# Returns
- `stat::Bool`: `false` if the peak is inside; `true` if not.
- `kpi::Float64`: the peak location x.
- `lpi::Float64`: the peak location y.
- `scrp::Float64`: the peak score.
"""
function find_score_peak(o::VTT, scr::Matrix{Float32}, kw::Int, lw::Int)
    # find the max score and its index (the last one, in case of ties, to match
    # the previous `findlast`-based behavior); NaN entries (which can occur when
    # `chk_mask` leaves too few visible samples in a window) are ignored.
    scrp = Float32(-Inf)
    lpi, kpi = 0, 0
    for k in axes(scr, 2), l in axes(scr, 1)
        v = @inbounds scr[l, k]
        if !isnan(v) && v >= scrp
            scrp = v
            lpi, kpi = l, k
        end
    end
    if kpi == 0
        return true, nothing, nothing, nothing
    end

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
- `o::VTT`: The tracking object.
- `tid::Array{Integer,Any}`: Tracking initial time indices.
- `x::Array{Float64,Any}`: Tracking initial template-center x location (index-based; non-integer for subgrid).
- `y::Array{Float64,Any}`: Tracking initial template-center y location (index-based; non-integer for subgrid).
- `vxg::Array{Float64,Any}=nothing`: First guess of vx (to search around it). Can be 0.
- `vyg::Array{Float64,Any}=nothing`: First guess of vy (to search around it). Can be 0.
- `out_subimage::Bool=false`: Whether output subimages.
- `out_score_ary::Bool=false`: Whether output score arrays.
- `to_missing::Bool=true`: Whether output missing values as `missing`.

# Returns
- `count::Vector{Integer}`: [len] The number of valid trajectory points for each initial
    template, including the initial one (so it ranges 0:ntrac+1; e.g. `ntrac+1` means every
    step succeeded).
- `status::Vector{Integer}`: [len] The reason tracking stopped for each initial template
    (`0` if it completed `ntrac` steps without stopping early). One of:
    - `0`: completed successfully (or not yet stopped)
    - `1`/`5`: the start/end `tid` for a step was out of range
    - `2`: the template subimage could not be read (out of bounds, or missing/masked data)
    - `3`: the template subimage failed the `Cth` contrast check
    - `4`: the template subimage failed the `chk_peak_inside` check
    - `6`: the sliding score could not be computed (out of bounds, or missing/masked data)
    - `7`: no valid score peak, or the peak was on the edge of the search window
    - `8`: the peak score was below `Sth0` (1st step) or `Sth1` (subsequent steps)
    - `9`: the velocity change along the trajectory exceeded `vxch`/`vych`
- `tid::Matrix{Float64}`: [ntrac+1, len] time index of the trajectories (tid0 and subsequent ones).
- `x::Matrix{Float64}`: [ntrac+1, len] x locations of the trajectories (x0 and derived ones).
- `y::Matrix{Float64}`: [ntrac+1, len] y locations of trajectories (x0 and derived ones).
- `vx::Matrix{Float64}`: [ntrac, len] Derived x-velocity.
- `vy::Matrix{Float64}`: [ntrac, len] Derived y-velocity.
- `score::Matrix{Float64}`: [ntrac, len] Scores along the trajectory (max values, possibly at subgrid).
- `zss::Array{Float32,4}`: [ntrac+1, nsy, nsx, len] optional, if non-`nothing`
    (Diagnosis output if wanted) The subimages along the track.
- `score_ary::Array{Float64,4}`: [ntrac, 2*iyhw+1, 2*ixhw+1, len] optional, if non-`nothing`
    (Diagnosis output if wanted) The entire scores.
"""
function trac(o::VTT, tid, x, y; vxg=nothing, vyg=nothing, out_subimage::Bool=false, out_score_ary::Bool=false, to_missing::Bool=true)
    !o.setuped && throw(ArgumentError("Need to call #setup in advance"))
    sh = size(x)
    if tid isa Integer
        tid = fill(tid, sh...)
    else
        size(tid) != sh && throw(ArgumentError("Shape miss-match (x)"))
    end

    size(y) != sh && throw(ArgumentError("Shape miss-match (y)"))

    if isnothing(vxg)
        vxg = fill(0.0, sh...)
    else
        size(vxg) != sh && throw(ArgumentError("Shape miss-match (vxg)"))
    end

    if isnothing(vyg)
        vyg = fill(0.0, sh...)
    else
        size(vyg) != sh && throw(ArgumentError("Shape miss-match (vyg)"))
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
        if !isnothing(zss)
            zss = Array{Union{Missing, Float32},4}(zss)
            zss[zss.==o.zmiss] .= missing
        end
        if !isnothing(score_ary)
            score_ary = Array{Union{Missing, Float32},4}(score_ary)
            score_ary[score_ary.==Float32(fmiss)] .= missing
        end
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


# Status codes for `trac`/`do_tracking`'s `status` output (`0` means the step succeeded):
const STATUS_TID_START_OUT_OF_RANGE = 1    # `tid` of the tracking start time is out of range
const STATUS_TEMPLATE_READ_FAILED = 2      # template subimage could not be read (out of bounds or missing data)
const STATUS_LOW_CONTRAST = 3              # template subimage failed the `Cth` contrast check
const STATUS_PEAK_NOT_INSIDE_TEMPLATE = 4  # template subimage failed the `chk_peak_inside` check
const STATUS_TID_END_OUT_OF_RANGE = 5      # `tid` of the tracking end time is out of range
const STATUS_SCORE_COMPUTATION_FAILED = 6  # sliding score could not be computed (out of bounds or missing data)
const STATUS_PEAK_NOT_FOUND = 7            # no valid score peak, or the peak was on the search-window edge
const STATUS_SCORE_BELOW_THRESHOLD = 8     # peak score below `Sth0` (1st step) or `Sth1` (subsequent steps)
const STATUS_VELOCITY_CHANGE_TOO_LARGE = 9 # velocity change along the trajectory exceeded `vxch`/`vych`

"""
    do_tracking(o, tid0, x0, y0, vx0, vy0, out_subimage, out_score_ary)

Conduct tracking (core).

# Arguments
- `o::VTT`: The tracking object.
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

    visible = trues(o.nsy, o.nsx) # placeholder; overwritten below when `chk_mask`.
    # Declared outside the `j` loop (like `zs0`) so it persists across steps
    # when `use_init_temp` is true and the template/visible is not recomputed.

    status = zeros(Int, len)

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
                status[m] = STATUS_TID_START_OUT_OF_RANGE
                continue
            end
            
            if j == 1 || !o.use_init_temp 
                if o.chk_mask
                    stat, zs0, visible = get_zsub_visible_subgrid(o, tidf, xcur, ycur)
                    if !stat
                        stat = sum(visible) < o.min_samples
                    end
                else
                    stat, zs0 = get_zsub_subgrid(o, tidf, xcur, ycur)
                end
            end

            if stat
                status[m] = STATUS_TEMPLATE_READ_FAILED
                continue
            end

            if o.chk_Cth
                if maximum(zs0) - minimum(zs0) < o.Cth
                    status[m] = STATUS_LOW_CONTRAST
                    continue
                end
            end

            if o.chk_peak_inside
                stat = chk_zsub_peak_inside(o, zs0)
                if stat
                    status[m] = STATUS_PEAK_NOT_INSIDE_TEMPLATE
                    continue
                end
            end

            if out_subimage
                if j == 1 || !o.use_init_temp
                    zss[j,:,:,m] = zs0
                else
                    stat, zsw = get_zsub_subgrid(o, tidf, xcur, ycur)
                    if stat
                        status[m] = STATUS_TEMPLATE_READ_FAILED
                        continue
                    end
                    zss[j,:,:,m] = zsw
                end
            end

            # inspect the tracking end time
            stat = inspect_t_index(o, tidl)
            if stat
                status[m] = STATUS_TID_END_OUT_OF_RANGE
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
                status[m] = STATUS_SCORE_COMPUTATION_FAILED
                continue
            end

            if out_score_ary
                score_ary[j,:,:,m] .= scr
            end

            stat, xp, yp, sp = find_score_peak(o, scr, kw, lw)
            if stat
                status[m] = STATUS_PEAK_NOT_FOUND
                continue
            end
            if ((j==1 && sp<o.Sth0) || (j>1 && sp<o.Sth1))
                status[m] = STATUS_SCORE_BELOW_THRESHOLD
                continue
            end
            xw = xp + kc - 1 - ixhw # next position (x-axis)
            yw = yp + lc - 1 - iyhw # next position (y-axis)
            vxw = (xw - x[j,m])/dt # velocity (x-axis)
            vyw = (yw - y[j,m])/dt # velocity (y-axis)
            if chk_vchange && j > 1
                stat = (abs(vxw - vx[j-1,m]) > o.vxch || abs(vyw - vy[j-1,m]) > o.vych)
                if stat
                    if j == 2 # no consecutive consistent result: invalidate the j==1 result too
                        count[m] = 0
                        tid[2,m] = imiss
                        x[2,m] = y[2,m] = o.fmiss
                        vx[1,m] = vy[1,m] = o.fmiss
                        score[1,m] = Float32(o.fmiss)
                    end
                    status[m] = STATUS_VELOCITY_CHANGE_TOO_LARGE
                    continue
                end
            end
            count[m] = j + 1 # of valid tracking
            tid[j+1,m] = tidl
            x[j+1,m] = xw
            y[j+1,m] = yw
            vx[j,m] = vxw
            vy[j,m] = vyw
            score[j,m] = Float32(sp)
            if out_subimage && j == o.ntrac # last sub image
                xcur = @inbounds x[j+1,m]
                ycur = @inbounds y[j+1,m]
                stat, zs0 = get_zsub_subgrid(o, tidl, xcur, ycur)
                if stat
                    status[m] = STATUS_TEMPLATE_READ_FAILED
                    continue
                end
                zss[j+1,:,:,m] .= zs0
            end
        end
    end   
    return count, status, tid, x, y, vx, vy, score, zss, score_ary
end
end
