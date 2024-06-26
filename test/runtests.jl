using VTTrac
using Test

using Statistics

@testset "VTTrac.jl" begin
    @testset "VTTrack_original" begin
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
        @test vtt.t == tax
        @test vtt.chk_zmiss == false


        @test_throws ArgumentError VTTrac.VTT(z, t=tax[begin:end-1])


        vtt = VTTrac.VTT(z, t=tax, zmiss=-999.0)
        @test vtt.nt == nt
        @test vtt.ny == ny
        @test vtt.nx == nx
        @test vtt.chk_zmiss == true
        @test vtt.setuped == false


        ntrac = nt-1
        nsx = 5
        nsy = 5
        @test_throws ArgumentError VTTrac.setup(vtt, nsx, nsy)
        @test_throws ArgumentError VTTrac.setup(vtt, nsx, nsy; vxhw=1.8, vyhw=1.8, ixhw=3, iyhw=3)
        

        VTTrac.setup(vtt, nsx, nsy; ixhw=3, iyhw=3, itstep=5)
        @test vtt.dtmean == 5.0

        VTTrac.setup(vtt, nsx, nsy; ixhw=3, iyhw=3)
        @test vtt.dtmean == 1

        @test vtt.vxhw == 2.0
        @test vtt.vyhw == 2.0

        VTTrac.setup(vtt, nsx, nsy; vxhw=1.8, vyhw=1.8, ntrac=ntrac, subgrid=false, subgrid_gaus=true, use_init_temp=false, score_method="xcor")
        @test vtt.ixhw == 3
        @test vtt.iyhw == 3
        @test vtt.peak_inside_th == -1.0f0
        @test vtt.Cth == -1.0f0


        n = 6
        tid0 = Vector{Int}(ones(n))
        x0 = 1 .+ [0:n-1;]*2.5 .+ 7.5
        y0 = 1 .+ [0:n-1;]*1.0 .+ 10.5
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test vtt.z == z # Check to see if the view is being written.
        @test count == fill(nt, n)
        @test tid == repeat([1:nt;]', 6)'
        @test mean(vx) == 1.0
        @test mean(vy) == 1.0
        @test size(count) == (n,)
        @test size(status) == (n,)
        @test size(tid) == (ntrac+1, n)
        @test size(x) == (ntrac+1, n)
        @test size(y) == (ntrac+1, n)
        @test size(vx) == (ntrac, n)
        @test size(vy) == (ntrac, n)
        @test size(score) == (ntrac, n)
        @test size(zss) == (ntrac+1, nsy, nsx, n)
        @test size(score_ary) == (ntrac, 2vtt.iyhw+1, 2vtt.ixhw+1, n)


        vtt.subgrid = true
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test 1.19 < mean(vx) < 1.21
        @test 1.19 < mean(vy) < 1.21


        vtt.chk_Cth = true
        vtt.Cth = 1.6
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test count !== fill(nt, n)
        @test 3 in status
        vtt.chk_Cth = false
        vtt.Cth = -1.0


        vtt.score_method = "ncov"
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test 1.19 < mean(vx) < 1.21
        @test 1.19 < mean(vy) < 1.21
        vtt.score_method = "xcor"


        vtt.chk_peak_inside = true
        vtt.peak_inside_th = 0.03
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test count !== fill(nt, n)
        vtt.chk_peak_inside = false
        vtt.peak_inside_th = -1.0


        vtt.use_init_temp = true
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test 1.19 < mean(vx) < 1.21
        @test 1.19 < mean(vy) < 1.21
        vtt.use_init_temp = false


        n1 = 2
        n2 = 3
        x0 = reshape(x0, n1, n2)
        y0 = reshape(y0, n1, n2)
        tid0 = reshape(tid0, n1, n2)
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test size(count) == (n1, n2)
        @test size(status) == (n1, n2)
        @test size(tid) == (ntrac+1, n1, n2)
        @test size(x) == (ntrac+1, n1, n2)
        @test size(y) == (ntrac+1, n1, n2)
        @test size(vx) == (ntrac, n1, n2)
        @test size(vy) == (ntrac, n1, n2)
        @test size(score) == (ntrac, n1, n2)
        @test size(zss) == (ntrac+1, nsy, nsx, n1, n2)
        @test size(score_ary) == (ntrac, 2vtt.iyhw+1, 2vtt.ixhw+1, n1, n2)
    end

    @testset "VTTrack_with_mask" begin
        nt = 20
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
        mask = z .>= maximum(z) # almost not masked
        vtt = VTTrac.VTT(z)
        @test vtt.chk_mask == false

        vtt = VTTrac.VTT(z, t=tax, mask=mask)
        @test vtt.visible == .!mask
        @test vtt.chk_mask == true

        ntrac = nt-1
        nsx = 5
        nsy = 5
        VTTrac.setup(vtt, nsx, nsy; vxhw=1.8, vyhw=1.8, ntrac=ntrac, subgrid=true, subgrid_gaus=true, use_init_temp=false, score_method="xcor")

        n = 6
        tid0 = Vector{Int}(ones(n))
        x0 = 1 .+ [0:n-1;]*2.5 .+ 7.5
        y0 = 1 .+ [0:n-1;]*1.0 .+ 10.5
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test vtt.visible == .!mask # Check to see if the view is being written.
        
        vtt.min_samples = vtt.nsx * vtt.nsy
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test any(ismissing.(score_ary))
        vtt.min_samples = 1

        mask = trues(size(z)) # all masked
        vtt.visible = .!(mask)
        vtt.chk_mask = true
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test any(ismissing.(score_ary))
        
        mask = falses(size(z)) # all visible
        vtt.visible = .!(mask)
        vtt.chk_mask = true
        count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(vtt, tid0, x0, y0, out_subimage=true, out_score_ary=true)
        @test !any(ismissing.(x))

    end
end
