using VTTrac

nt = 10
ny = 100
nx = 100
tax = Vector{Float64}([0:nt-1;])
xax = [0:nx-1;]
yax = [0:ny-1;]
xg = permutedims(repeat(xax', outer=(length(yax),1,length(tax))), (3,1,2))
yg = permutedims(repeat(yax, outer=(1,length(xax),length(tax))), (3,1,2))
tg = repeat(tax, outer=(1, length(yax), length(xax)))
k = 2pi/10
cx = 1.2
cy = 1.2
z = sin.(k*(xg-cx*tg)) .* cos.(k*(yg-cy*tg))
z = Array{Float32,3}(z)
vtt = VTTrac.VTT(z)

ntrac = nt-1
nsx = 5
nsy = 5

VTTrac.setup(vtt, nsx, nsy; vxhw=1.8, vyhw=1.8, ntrac=ntrac, subgrid=false, subgrid_gaus=true, use_init_temp=false, score_method="xcor")

n = 6
tid0 = Vector{Int}(ones(n))
x0 = 1 .+ [0:n-1;]*2.5 .+ 7.5
y0 = 1 .+ [0:n-1;]*1.0 .+ 10.5
count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
