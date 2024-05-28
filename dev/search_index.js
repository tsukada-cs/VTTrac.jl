var documenterSearchIndex = {"docs":
[{"location":"references/","page":"References","title":"References","text":"CurrentModule = VTTrac","category":"page"},{"location":"references/","page":"References","title":"References","text":"setup\nset_ixyhw_from_v\nset_ixyhw_directly\ntrac","category":"page"},{"location":"references/#VTTrac.setup","page":"References","title":"VTTrac.setup","text":"setup(o, nsx, nsy, vxhw, vyhw, [ixhw, iyhw, subgrid, subgrid_gaus,\n    itstep, ntrac, score_method, Sth0, Sth1, vxch, vych,\n    peak_inside, peak_inside_th, Cth, use_init_temp, min_samples])\n\nSetup for tracking.\n\nArguments\n\no::VTT: The object.\nnsx::Integer, nsy::Integer: Submimage x & y sizes (x:1st, y:2nd dim).\nvxch::Union{Real, Nothing}, vyhw::Union{Real, Nothing}: (either v[xy]hw or i[xy]hw are MANDATORY).   the dimensions along which to perform the computation.   search velocity range half sizes to set i[xy]hw.   Seach at least to cover +-v?hw around the first guess or previous step.   (the result can be outside the range.)\nixhw::Union{Int, Nothing}, iyhw::Union{Int, Nothing}: (either v[xy]hw or i[xy]hw are MANDATORY)   Max displacement fro template match (can be set indirecly through v[xy]hw).\nsubgrid::Bool=true: Whether to conduct subgrid tracking.\nsubgrid_gaus::Bool=true: Whether subgrid peak finding is by gaussian.\nitstep::Integer=1: Step of t's used (skip if >1).\nntrack::Integer=2: Max tracking times from initial loc.\nscore_method::String=\"xcor\": \"xcor\" for cross-correlation, \"ncov\" for normalized covariance.\nSth0::AbstractFloat=0.8: The minimum score required for the 1st tracking.\nSth1::AbstractFloat=0.7: The minimum score required for subsequent tracking.\nvxch::Union{Real, Nothing}=nothing: If non-nothing, the max tolerant vx   change between two consecutive tracking.\nvych::Union{Real, Nothing}=nothing: If non-nothing, the max tolerant vy   change between two consecutive tracking.\npeak_inside_th::Union{Real, Nothing}=nothing: If non-nothing, an initial template is used only when   it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.\nCth::Union{Real, Nothing}=nothing: If non-nothing, an initial template is used only when    it has a difference in max and min greater than its value.\nmin_samples::Int=1: Minimum number of visible values to calculate score when chk_mask is true.\n\n\n\n\n\n","category":"function"},{"location":"references/#VTTrac.set_ixyhw_from_v","page":"References","title":"VTTrac.set_ixyhw_from_v","text":"set_ixyhw_from_v(o, vxch, vyxh)\n\nSets the tracking parameters i[xy]hw from velocities (v[xy]hh).\n\nArguments\n\no::VTT: The object.\nvxhw::Float64: The range over which vx is searched around initial guess.\nvyhw::Float64: The range over which vy is searched around initial guess.\n\n```\n\n\n\n\n\n","category":"function"},{"location":"references/#VTTrac.set_ixyhw_directly","page":"References","title":"VTTrac.set_ixyhw_directly","text":"set_ixyhw_from_v(o, ixch, iyxh)\n\nSets the tracking parameters i[xy]hw.\n\nArguments\n\no::VTT: The object.\nixhw::Float64: The range over which next x is searched around initial guess.\niyhw::Float64: The range over which next y is searched around initial guess.\n\n\n\n\n\n","category":"function"},{"location":"references/#VTTrac.trac","page":"References","title":"VTTrac.trac","text":"trac(o, tid, x, y[, vxg, vyg, out_subimage, out_score_ary])\n\nConduct tracking.\n\nArguments\n\no::VTT: The traking oject.\ntid::Array{Integer,Any}: Tracking initial time indices.\nx::Array{Float64,Any}: Tracking initial template-center x location (index-based; non-integer for subgrid).\ny::Array{Float64,Any}: Tracking initial template-center y location (index-based; non-integer for subgrid).\nvxg::Array{Float64,Any}=nothing: First guess of vx (to search around it). Can be 0.\nvyg::Array{Float64,Any}=nothing: First guess of vy (to search around it). Can be 0.\nout_subimage::Bool=false: Whether output subimages.\nout_score_ary::Bool=false: Whether output score arrays.\nto_missing::Bool=true: Whether output missing values as missing.\n\nReturns\n\ncount::Vector{Integer}: [len] The number of successful tracking for each initial template.\ntid::Matrix{Float64}: [ntrac+1, len] time index of the trajectories (tid0 and subsequent ones).\nx::Matrix{Float64}: [ntrac+1, len] x locations of the trajectories (x0 and derived ones).\ny::Matrix{Float64}: [ntrac+1, len] y locations of trajectories (x0 and derived ones).\nvx::Matrix{Float64}: [ntrac, len] Derived x-velocity.\nvy::Matrix{Float64}: [ntrac, len] Derived y-velocity.\nscore::Matrix{Float64}: [ntrac, len] Scores along the trajectory (max values, possibly at subgrid).\nzss::Array{Float32,4}: [nsx, nsy, ntrac+1, len] optional, if non-nothing   (Diagnosis output if wanted) The subimages along the track.\nscore_arry::Array{Float64,4}: [(x-sliding size, y-sliding size, ntrac+1, len] optional, if non-nothing   (Diagnosis output if wanted) The entire scores.\n\n\n\n\n\n","category":"function"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"CurrentModule = VTTrac","category":"page"},{"location":"howtouse/#How-to-Use","page":"How to Use","title":"How to Use","text":"","category":"section"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[1] Import VTTrac","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    using VTTrac","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[2] Make a sample time series data","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    nt = 10\n    ny = 100\n    nx = 100\n    tax = Vector{Float64}([0:nt-1;])\n    xax = [0:nx-1;]\n    yax = [0:ny-1;]\n    xg = permutedims(repeat(xax', outer=(length(yax),1,length(tax))), (3,1,2))\n    yg = permutedims(repeat(yax, outer=(1,length(xax),length(tax))), (3,1,2))\n    tg = repeat(tax, outer=(1, length(yax), length(xax)))\n    k = 2pi/10\n    cx = 1.2\n    cy = 1.2\n    z = sin.(k*(xg-cx*tg)) .* cos.(k*(yg-cy*tg))\n    z = Array{Float32,3}(z)","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[3] Initialize VTTrac.VTT object","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    vtt = VTTrac.VTT(z)","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[4] Setup tracking","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    ntrac = nt-1 # number of tracking\n    nsx = 5 # template size in x axis\n    nsy = 5 # template size in y axis\n\n    VTTrac.setup(vtt, nsx, nsy; vxhw=1.8, vyhw=1.8, ntrac=ntrac, subgrid=false, subgrid_gaus=true, use_init_temp=false, score_method=\"xcor\")","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[5] Set initial template positions","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    n = 6\n    tid0 = Vector{Int}(ones(n))\n    x0 = 1 .+ [0:n-1;]*2.5 .+ 7.5\n    y0 = 1 .+ [0:n-1;]*1.0 .+ 10.5","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"[6] Conduct tracking","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"    count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)","category":"page"},{"location":"howtouse/","page":"How to Use","title":"How to Use","text":"To see details of outputs, see VTTrac.trac.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VTTrac","category":"page"},{"location":"#VTTrac","page":"Home","title":"VTTrac","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for VTTrac.","category":"page"},{"location":"#VTTrac:-Velocimetry-by-Template-Tracking","page":"Home","title":"VTTrac: Velocimetry by Template Tracking","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This library provides the julia-language implementation for VTTrac (https://github.com/thorinouchi/VTTrac). It does not use module variables, so it should be good for parallel execution.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm used in this library is the simple template matching of PIV (particle image velocimetry) for monochromatic image-like data, but the matching is conducted multiple times in a Lagrangian manner as in PTV (particle tracking velocimetry) over a number of times specified by the parameter named ntrac. The default scoring method for template matching is the cross correlation coefficient, as in the basic PIV. Both forward and backward tracking is available. Use the parameter itstep; tracking is backward along time sequence, if it is negative.","category":"page"},{"location":"#Available-scoring-methods","page":"Home","title":"Available scoring methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"XCOR: cross-correlation, cov(x',y')/sig(x)/sig(y), where x represents the template sub-image and y represents the target sub-image that are slided.\nNCOV: normalized covariance, cov(x',y')/sig(x)^2: covariance normalized the variance of the template sub-image x.","category":"page"},{"location":"#Screening","page":"Home","title":"Screening","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A check (result screening) based on velocity change along trajectory is available (by using the threshold parameter named vxch and vych), so it is recommended to always set ntrac >= 2. Further screening is available for initial templates (e.g., in terms of the complexity and contrast) and the quality of the results (e.g., score threshold); see the source code.","category":"page"},{"location":"#Dimensions","page":"Home","title":"Dimensions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Spatial coordinates are based on array indices, with the distance between adjacent grid points always being 1, so they are non-dimensional. The velocities are based on non-dimensional spacial displacement over time difference, where, time can either be dimensional or non-dimensional.","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"VTTrac by Takeshi Horinouchi: https://github.com/thorinouchi/VTTrac\npyVTTrac by Taiga Tsukada: https://github.com/tsukada-cs/pyVTTrac","category":"page"}]
}
