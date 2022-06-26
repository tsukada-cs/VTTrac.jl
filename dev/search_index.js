var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = VTTrac","category":"page"},{"location":"#VTTrac","page":"Home","title":"VTTrac","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VTTrac.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VTTrac]","category":"page"},{"location":"#VTTrac.chk_zmiss_region-Tuple{VTT, Any, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.chk_zmiss_region","text":"chk_zmiss_region(o, tid, k0, k1, l0, l1)\n\nCheck if there is data missing in the specified region at tid.\n\nReturns\n\nstat::Bool: false if there is no data missing, true if not.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.chk_zsub_peak_inside-Tuple{VTT, Any}","page":"Home","title":"VTTrac.chk_zsub_peak_inside","text":"chk_zsub_peak_inside(o, zs)\n\nCheck whether the template subimage is peaked (maximized or minimized) inside, and the peak is conspicuous enough, having a difference from the  max or min on the sides greater than peakinsideth*(insidemax - insidemin).\n\nCaution\n\nIf o.peak_inside_th < 0, no checking is conducted.\n\nReturns\n\nstat::Bool: false if passed the check, true if not.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.do_tracking-Tuple{VTT, Any, Any, Any, Any, Any, Bool, Bool}","page":"Home","title":"VTTrac.do_tracking","text":"do_tracking(o, tid0, x0, y0, vx0, vy0, out_subimage, out_score_ary)\n\nConduct tracking (core).\n\nArguments\n\no::VTT: The traking oject.\ntid0::Vector{Integer}: Tracking initial time indices.\nx0::Vector{Float64}: Tracking initial template-center x location (index-based; non-integer for subgrid).\ny0::Vector{Float64}: Tracking initial template-center y location (index-based; non-integer for subgrid).\nvx0g::Vector{Float64}: First guess of vx (to search around it). Can be 0.\nvy0g::Vector{Float64}: First guess of vy (to search around it). Can be 0.\nout_subimage::Bool: Whether output subimages.\nout_score_ary::Bool: Whether output score arrays.\n\nReturns\n\ncount::Array{Float64,2}: (len: len) The number of successful tracking for each initial template.\ntid::Array{Float64,2}: (len: (ntrac+1)*len) Time index of the trajectories (tid0 and subsequent ones).\nx::Array{Float64,2}: (len: (ntrac+1)*len) x locations of the trajectories (x0 and derived ones).\ny::Array{Float64,2}: (len: (ntrac+1)*len) y locations of trajectories (x0 and derived ones).\nvx::Array{Float64,2}: (len: ntrac*len) Derived x-velocity.\nvy::Array{Float64,2}: (len: ntrac*len) Derived y-velocity.\nscore::Array{Float64,2}: (len: ntrac*len)  Scores along the trajectory (max values, possibly at subgrid).\nzss::Array{Float64,4}: (optional, if non-nothing) (Diagnosis output if wanted) The subimages along the track (1D pointer for 4D array; nsx * nsy * (ntrac+1) * len.\nscore_arry::Array{Float64,4}: (optional, if non-nothing) (Diagnosis output if wanted) the entire scores (1D pointer for 4D array; (x-sliding size) * (y-sliding size) * (ntrac+1) * len.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.find_score_peak-Tuple{VTT, Any, Int64, Int64}","page":"Home","title":"VTTrac.find_score_peak","text":"find_score_peak(o, scr, kw, lw)\n\nFind the score peak and its location.\n\nReturns\n\nstat::Bool: false if the peak is inside; true if not.\nkpi::Float64: the peak location x.\nlpi::Float64: the peak location y.\nscrp::Float64: the peak score.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.find_subgrid_peak_5pt_epara-NTuple{5, Real}","page":"Home","title":"VTTrac.find_subgrid_peak_5pt_epara","text":"find_subgrid_peak_5pt_epara(c, l, r, b, t)\n\nFind subgrid peak from 5 points with elliptic paraboloid.\n\nInput c(enter), l(eft), r(ight), b(ottom), t(op) at (0,0), (-1,0), (0,1), (0,-1), (0,1), respectively. c must be greater than any of l,r,b,t.\n\nEquation: z = -p(x-x0)^2 + -q(y-y0)^2 + r             = -p x^2 + 2p x0 x - q y^2 + 2q y0 y + c\n\n, where u = p x0, v = q y0, c = -p x0^2 - q y0^2 + r\n\nReturns\n\nstat::Bool: false if find peak successfully, true if not.\nx0::float: x of peak the location.\ny0::float: y of peak the location.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.find_subgrid_peak_5pt_gaus-NTuple{5, Real}","page":"Home","title":"VTTrac.find_subgrid_peak_5pt_gaus","text":"find_subgrid_peak_5pt_gaus(c, l, r, b, t)\n\nFind subgrid peak from 5 points by interpolating with a 2D gaussian (for positive scores).\n\nIt is simply a log-version of the elliptic-paraboloid method (find_subgrid_peak_5pt_epara). It appears that this method is preferred in many PIVs over the elliptic-paraboloid method (find_subgrid_peak_5pt_epara). In a test, there were not much difference between them, though.\n\nInput c(enter), l(eft), r(ight), b(ottom), t(op) at (0,0), (-1,0), (0,1), (0,-1), (0,1), respectively. all of them must be positive. c must be greater than any of l,r,b,t.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_score-Tuple{VTT, Any, Any, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.get_score","text":"Conduct template matching driver.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_score_ncov-Tuple{VTT, Any, Int64, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.get_score_ncov","text":"get_score_ncov(o, x, tid, k0, k1, l0, l1)\n\nConduct template matching, scoring by normalized covariance.\n\nReturns\n\nstat::Bool: false if passed the check, true if not.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_score_xcor-Tuple{VTT, Any, Int64, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.get_score_xcor","text":"get_score_xcor(o, x, tid, k0, k1, l0, l1)\n\nConduct template matching, scoring by cross-correlation.\n\nReturns\n\nstat::Bool: false if passed the check, true if not.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_zsub-Tuple{VTT, Int64, Int64, Int64}","page":"Home","title":"VTTrac.get_zsub","text":"get_zsub(o, tid, xi, yi)\n\nRead out a template subimage from the image at tid.\n\nThe sub-image positons are specified at its center. (If the sub-image size is even, with one more pix on the \"left\" / \"bottom\", by starting from the index xi-nsx/2, yi-nsy/2).\n\nReturns\n\nstats::Bool: false if successful (specified region is valid and, if chk_zmiss, \n\nno data missing), true if not.\n\nzsub::Array{Float32,2}: Subimage at (x,y) = (xi, yi). \n\nSee Also\n\nget_zsub_view\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_zsub_subgrid-Tuple{VTT, Int64, Float64, Float64}","page":"Home","title":"VTTrac.get_zsub_subgrid","text":"get_zsub_subgrid(o, tid, x, y)\n\nRead out a template submimage from the image at tid. Possibly at subgrid: Linearly interpolated, if x or y has deviation from integer (bilinear if x and y). Efficient: no unnecessary read-out is made.\n\nReturns\n\nstats::Bool: false if successful (specified region is valid and, if chk_zmiss, \n\nno data missing), true if not.\n\nzsubg::Array{Float32,2}: Subimage.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.get_zsub_view-Tuple{VTT, Int64, Int64, Int64}","page":"Home","title":"VTTrac.get_zsub_view","text":"get_zsub_view(o, tid, xi, yi)\n\nLike get_zsub, but returns a view.\n\nNotes\n\nThis function returns a view of subimage. get_zsub returns a copy of subimage.\n\nSee Also\n\nget_zsub\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.inspect_t_index-Tuple{VTT, Int64}","page":"Home","title":"VTTrac.inspect_t_index","text":"To check whether a time index is valid. Returns false if valid, true if not.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.print_ary2d-Tuple{Any}","page":"Home","title":"VTTrac.print_ary2d","text":"print a 2D double array\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.roundInt-Tuple{Any}","page":"Home","title":"VTTrac.roundInt","text":"round to Int like C/C++\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.set_basic!-Tuple{VTT, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.set_basic!","text":"set_basic!(o, nsx, nsy, itstep, ntrac)\n\nSets basic (mandatory) tracking parameters. Also sets default vals for optional params.\n\nArguments\n\no::VTT: The object.\nnsx::Integer: The template subimage size (x).\nnsy::Integer: The template subimage size (y).\nitstep::Integer: Index-based time step for tracking (can be negative).\nntrac::Integer: Number of times for each initial template is tracked.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.set_ixyhw_directly!-Tuple{VTT, Int64, Int64}","page":"Home","title":"VTTrac.set_ixyhw_directly!","text":"set_ixyhw_from_v!(o, ixch, iyxh)\n\nSets the tracking parameters i[xy]hw.\n\nArguments\n\no::VTT: The object.\nixhw::Float64: The range over which next x is searched around initial guess.\niyhw::Float64: The range over which next y is searched around initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.set_ixyhw_from_v!-Tuple{VTT, Float64, Float64}","page":"Home","title":"VTTrac.set_ixyhw_from_v!","text":"set_ixyhw_from_v!(o, vxch, vyxh)\n\nSets the tracking parameters i[xy]hw from velocities (v[xy]hh).\n\nArguments\n\no::VTT: The object.\nvxhw::Float64: The range over which vx is searched around initial guess.\nvyhw::Float64: The range over which vy is searched around initial guess.\n\n```\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.set_optional!-Tuple{VTT, Bool, Bool, String, AbstractFloat, AbstractFloat, Real, Real, Float64, Float64, Bool}","page":"Home","title":"VTTrac.set_optional!","text":"set_optional!(o, subgrid, subgrid_gaus, score_method, score_th0, score_th1, peak_inside_th, min_contrast, vxch, vych, use_init_temp)\n\nSets optional tracking parameters.\n\nArguments\n\no::VTT: The object.\nsubgrid::Bool: Whether to conduct subgrid tracking.\nsubgrid_gaus::Bool: Whether subgrid peak finding is by gaussian.\nscore_method::String: Scoring method (such as xcor for cross-correlation).\nscore_th0::AbstractFloat: (Result screening parameter) Minimum score required for the first-time tracking.\nscore_th1::AbstractFloat: (Result screening parameter) Minimum score required for the subsequent tracking.\npeak_inside_th::Real: An initial template is used only when   it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.\nmin_contrast::Real: An initial template is used only when    it has a difference in max and min greater than its value.\nvxch::Float64: (Result screening parameter) If positive, tracking result is rejected if the vx    chnages along trajecty greather than this value (thus used only when ntrac>=2). As a special case,    if the result of the second tracking is rejected, the first one is also rejected, since there is    no consecutive consistent result in this case.\nvych::Float64: (Result screening parameter) As vxch but for the y-component.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.setup-Tuple{VTT, Int64, Int64}","page":"Home","title":"VTTrac.setup","text":"setup(o, nsx, nsy, vxhw, vyhw, [ixhw, iyhw, subgrid, subgrid_gaus,\n    itstep, ntrac, score_method, score_th0, score_th1, vxch, vych,\n    peak_inside, peak_inside_th, min_contrast, use_init_temp])\n\nSetup for tracking.\n\nArguments\n\no::VTT: The object.\nnsx::Integer, nsy::Integer: Submimage x & y sizes (x:1st, y:2nd dim).\nvxch::Union{Real, Nothing}, vyhw::Union{Real, Nothing}: (either v[xy]hw or i[xy]hw are MANDATORY).   the dimensions along which to perform the computation.   search velocity range half sizes to set i[xy]hw.   Seach at least to cover +-v?hw around the first guess or previous step.   (the result can be outside the range.)\nixhw::Union{Int, Nothing}, iyhw::Union{Int, Nothing}: (either v[xy]hw or i[xy]hw are MANDATORY)   Max displacement fro template match (can be set indirecly through v[xy]hw).\nsubgrid::Bool=true: Whether to conduct subgrid tracking.\nsubgrid_gaus::Bool=true: Whether subgrid peak finding is by gaussian.\nitstep::Integer=1: Step of t's used (skip if >1).\nntrack::Integer=2: Max tracking times from initial loc.\nscore_method::String=\"xcor\": \"xcor\" for cross-correlation, \"ncov\" for normalized covariance.\nscore_th0::AbstractFloat=0.8: The minimum score required for the 1st tracking.\nscore_th1::AbstractFloat=0.7: The minimum score required for subsequent tracking.\nvxch::Union{Real, Nothing}=nothing: If non-nothing, the max tolerant vx   change between two consecutive tracking.\nvych::Union{Real, Nothing}=nothing: If non-nothing, the max tolerant vy   change between two consecutive tracking.\npeak_inside_th::Union{Real, Nothing}=nothing: If non-nothing, an initial template is used only when   it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.\nmin_contrast::Union{Real, Nothing}=nothing: If non-nothing, an initial template is used only when    it has a difference in max and min greater than its value.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.sliding_ncov-Tuple{VTT, Any, Any, Int64, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.sliding_ncov","text":"sliding_ncov(o, sigx, xd, tid, k0, k1, l0, l1)\n\nSliding normalized covariance between the sugimage and image at tid.\n\nNormalization is done by the sigma of the fist image : cov(x',y')/sigx^2 (in contrast to cov(x',y')/sigx/sigy in the correlation coefficient).\n\nReturns\n\nstat::Bool: false if all the relevant data and regions are valid, so all the scores   (ncov) are defined at all tested center locations; true if not.\nscr::Array{Float64,2}: Score array.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.sliding_xcor-Tuple{VTT, Any, Any, Int64, Int64, Int64, Int64, Int64}","page":"Home","title":"VTTrac.sliding_xcor","text":"sliding_xcor(o, sigx, xd, tid, k0, k1, l0, l1)\n\nSliding cross-correlation between the sugimage and image at tid.\n\nReturns\n\nstat::Bool: false if all the relevant data and regions are valid, so all the scores   (xcor) are defined at all tested center locations; true if not.\nscr::Array{Float64,2}: Score array.\n\n\n\n\n\n","category":"method"},{"location":"#VTTrac.trac-Tuple{VTT, Any, Any, Any}","page":"Home","title":"VTTrac.trac","text":"trac(o, tid, x, y[, vxg, vyg, out_subimage, out_score_ary])\n\nConduct tracking.\n\nArguments\n\no::VTT: The traking oject.\ntid::Array{Integer,Any}: Tracking initial time indices.\nx::Array{Float64,Any}: Tracking initial template-center x location (index-based; non-integer for subgrid).\ny::Array{Float64,Any}: Tracking initial template-center y location (index-based; non-integer for subgrid).\nvxg::Array{Float64,Any}=nothing: First guess of vx (to search around it). Can be 0.\nvyg::Array{Float64,Any}=nothing: First guess of vy (to search around it). Can be 0.\nout_subimage::Bool=false: Whether output subimages.\nout_score_ary::Bool=false: Whether output score arrays.\nto_missing::Bool=true: Whether output missing values as missing.\n\nReturns\n\ncount::Vector{Integer}: [len] The number of successful tracking for each initial template.\ntid::Array{Float64,2}: [ntrac+1, len] time index of the trajectories (tid0 and subsequent ones).\nx::Array{Float64,2}: [ntrac+1, len] x locations of the trajectories (x0 and derived ones).\ny::Array{Float64,2}: [ntrac+1, len] y locations of trajectories (x0 and derived ones).\nvx::Array{Float64,2}: [ntrac, len] Derived x-velocity.\nvy::Array{Float64,2}: [ntrac, len] Derived y-velocity.\nscore::Array{Float64,2}: [ntrac, len] Scores along the trajectory (max values, possibly at subgrid).\nzss::Array{Float64,4}: nsx, nsy, ntrac+1, len   (Diagnosis output if wanted) The subimages along the track.\nscore_arry::Array{Float64,4}: (x-sliding size, y-sliding size, ntrac+1, len   (Diagnosis output if wanted) The entire scores.\n\n\n\n\n\n","category":"method"}]
}
