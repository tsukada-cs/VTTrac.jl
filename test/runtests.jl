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
        @test eltype(status) == Int
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
        
        # Test chk_zsub_peak_inside_function
        vtt.peak_inside_th = 0.1f0
        zs_peak = [
            0.0f0 0.0f0 0.0f0; 
            0.0f0 1.0f0 0.0f0; 
            0.0f0 0.0f0 0.0f0
        ]
        @test VTTrac.chk_zsub_peak_inside(vtt, zs_peak) == false
        zs_flat = ones(Float32, 3, 3)
        @test VTTrac.chk_zsub_peak_inside(vtt, zs_flat) == true
        zs_trough = [
            1.0f0 1.0f0 1.0f0; 
            1.0f0 0.0f0 1.0f0; 
            1.0f0 1.0f0 1.0f0
        ]
        @test VTTrac.chk_zsub_peak_inside(vtt, zs_trough) == false

        # to_missing=true with out_subimage=false or out_score_ary=false must not error
        x0_v = vec(x0); y0_v = vec(y0); tid0_v = vec(tid0)
        @test_nowarn VTTrac.trac(vtt, tid0_v, x0_v, y0_v, out_subimage=false, out_score_ary=false)
        @test_nowarn VTTrac.trac(vtt, tid0_v, x0_v, y0_v, out_subimage=true,  out_score_ary=false)
        @test_nowarn VTTrac.trac(vtt, tid0_v, x0_v, y0_v, out_subimage=false, out_score_ary=true)

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

    @testset "VTTrack_bugfixes" begin
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

        @testset "use_init_temp with mask does not error" begin
            mask = falses(size(z))
            mask[1,1,1] = true # ensure `chk_mask` is triggered (any(mask))
            vtt = VTTrac.VTT(z, t=tax, mask=mask)
            VTTrac.setup(vtt, 5, 5; vxhw=1.8, vyhw=1.8, ntrac=3, use_init_temp=true)
            n = 3
            tid0 = ones(Int, n)
            x0 = [20.0, 30.0, 40.0]
            y0 = [20.0, 30.0, 40.0]
            @test_nowarn VTTrac.trac(vtt, tid0, x0, y0)
        end

        @testset "final out_subimage is read from the correct (final) time" begin
            vtt = VTTrac.VTT(z, t=tax)
            VTTrac.setup(vtt, 5, 5; vxhw=1.8, vyhw=1.8, ntrac=3, subgrid=true)
            n = 4
            tid0 = ones(Int, n)
            x0 = [20.0, 30.0, 40.0, 50.0]
            y0 = [20.0, 25.0, 30.0, 35.0]
            count, status, tid, x, y, vx, vy, score, zss, _ = VTTrac.trac(vtt, tid0, x0, y0; out_subimage=true, to_missing=false)
            @test all(count .== size(tid, 1)) # sanity: every trajectory tracked fully
            for m in 1:n
                last_tid = tid[end, m]
                stat, zs_expected = VTTrac.get_zsub_subgrid(vtt, last_tid, x[end,m], y[end,m])
                @test !stat
                @test zss[end, :, :, m] == zs_expected
            end
        end

        @testset "vxch/vych rollback invalidates the derived point, not the initial one" begin
            nt3 = 3
            ny3 = 30
            nx3 = 100
            tax3 = Float64.(0:nt3-1)
            xax3 = 0:nx3-1
            yax3 = 0:ny3-1
            k3 = 2pi/10
            shiftx = [0.0, 1.0, 4.0] # cumulative x-shift per frame; speed jumps from 1 to 3
            z3 = Array{Float32,3}(undef, nt3, ny3, nx3)
            for i in 1:nt3, iy in 1:ny3, ix in 1:nx3
                z3[i,iy,ix] = sin(k3*(xax3[ix]-shiftx[i])) * cos(k3*yax3[iy])
            end
            vtt3 = VTTrac.VTT(z3, t=tax3)
            VTTrac.setup(vtt3, 5, 5; ixhw=6, iyhw=3, ntrac=2, subgrid=false, vxch=1.5, vych=1.5, Sth0=0.5, Sth1=0.5)
            count, status, tid, x, y, vx, vy, score, _, _ = VTTrac.trac(vtt3, [1], [50.0], [15.0]; to_missing=false)
            @test status[1] == 9
            @test count[1] == 0
            # the initial point is always retained, regardless of downstream failure
            @test x[1,1] == 50.0
            @test y[1,1] == 15.0
            @test tid[1,1] == 1
            # the now-untrusted j==1 result must be invalidated, not the initial point
            @test x[2,1] == vtt3.fmiss
            @test y[2,1] == vtt3.fmiss
            @test tid[2,1] == vtt3.imiss
            @test vx[1,1] == vtt3.fmiss
            @test vy[1,1] == vtt3.fmiss
            @test score[1,1] == Float32(vtt3.fmiss)
        end

        @testset "sliding_xcor is numerically stable for offset data" begin
            # A large mean offset relative to the pattern's variance is known to make a
            # naive incremental "E[y^2] - E[y]^2" variance update lose precision (and, in
            # the old code, corrupt all subsequent columns once it went negative). Compare
            # against a brute-force, per-window (numerically stable) reference.
            ny4, nx4 = 8, 40
            offset = 1000.0f0
            z4 = Array{Float32,3}(undef, 2, ny4, nx4)
            for iy in 1:ny4, ix in 1:nx4
                v = offset + 3.0f0*sin(Float32(2pi*(ix-1)/9))
                z4[1,iy,ix] = v
                z4[2,iy,ix] = v
            end
            vtt4 = VTTrac.VTT(z4, t=[0.0, 1.0])
            vtt4.nsx = 3; vtt4.nsy = 3; vtt4.chk_zmiss = false

            xtmpl = z4[1, 3:5, 3:5]
            xm = mean(xtmpl)
            xd = Float32.(xtmpl .- xm)
            sigx = Statistics.stdm(xtmpl, xm, corrected=false)

            k0, k1, l0, l1 = 5, 35, 3, 5
            stat, scr = VTTrac.sliding_xcor(vtt4, sigx, xd, 1, k0, k1, l0, l1)
            @test !stat

            nsx2, nsy2 = 1, 1
            nsxy = 9
            ref = zeros(Float32, l1-l0+1, k1-k0+1)
            for (li, l) in enumerate(l0:l1), (ki, k) in enumerate(k0:k1)
                sub = z4[1, l-nsy2:l-nsy2+2, k-nsx2:k-nsx2+2]
                ymean = mean(sub)
                yd = sub .- ymean
                vyy = sum(yd.^2)/nsxy
                vxy = sum(xd .* yd)/nsxy
                ref[li,ki] = vxy/sqrt(vyy)/sigx
            end
            @test maximum(abs.(scr .- ref)) < 1.0f-3
        end

        @testset "masked ncov uses the same cov(x,y)/var(x) normalization as unmasked ncov" begin
            ny5, nx5 = 40, 40
            nsx5, nsy5 = 5, 5
            z5 = Array{Float32,3}(undef, 2, ny5, nx5)
            for iy in 1:ny5, ix in 1:nx5
                v = 300.0f0 + Float32(sin(0.4*ix)*cos(0.3*iy))
                z5[1,iy,ix] = v
                z5[2,iy,ix] = v
            end
            vtt5 = VTTrac.VTT(z5)
            vtt5.nsx = nsx5; vtt5.nsy = nsy5; vtt5.chk_zmiss = false; vtt5.min_samples = 1
            nsx52, nsy52 = div(nsx5,2), div(nsy5,2)
            kc, lc = 20, 20
            ixhw, iyhw = 5, 5
            xtmpl = z5[1, lc-nsy52:lc-nsy52+nsy5-1, kc-nsx52:kc-nsx52+nsx5-1]
            xm = mean(xtmpl)
            xd = Float32.(xtmpl .- xm)
            sigx = Statistics.stdm(xtmpl, xm, corrected=false)
            k0, k1, l0, l1 = kc-ixhw, kc+ixhw, lc-iyhw, lc+iyhw

            # unmasked reference score
            vtt5.visible = trues(2, ny5, nx5)
            _, scr_unmasked = VTTrac.sliding_ncov(vtt5, sigx, xd, 2, k0, k1, l0, l1)

            # force the masked (loop) code path by making one pixel invisible far outside
            # the tested windows, while the windows themselves remain fully visible
            vis2 = trues(ny5, nx5)
            vis2[1,1] = false
            vtt5.visible = cat(trues(1, ny5, nx5), reshape(vis2, 1, ny5, nx5), dims=1)
            vis_tmpl = trues(nsy5, nsx5)
            _, scr_masked = VTTrac.get_score_ncov_with_visible(vtt5, xtmpl, vis_tmpl, 2, k0, k1, l0, l1)

            @test scr_unmasked == scr_masked
        end

        @testset "setup accepts Integer/Real broadly, not just Int/Float64" begin
            vtt = VTTrac.VTT(z, t=tax)
            @test_nowarn VTTrac.setup(vtt, 5, 5; vxhw=2, vyhw=2) # Int, not Float64
            vtt2 = VTTrac.VTT(z, t=tax)
            @test_nowarn VTTrac.setup(vtt2, 5, 5; ixhw=Int32(3), iyhw=Int32(3))
        end

        @testset "trac accepts non-Int64 scalar tid" begin
            vtt = VTTrac.VTT(z, t=tax)
            VTTrac.setup(vtt, 5, 5; ixhw=3, iyhw=3)
            @test_nowarn VTTrac.trac(vtt, Int32(1), [20.0], [20.0])
        end

        @testset "setup validates score_method" begin
            vtt = VTTrac.VTT(z, t=tax)
            @test_throws ArgumentError VTTrac.setup(vtt, 5, 5; ixhw=3, iyhw=3, score_method="bogus")
        end

        @testset "setup validates subimage size when peak_inside_th is set" begin
            vtt = VTTrac.VTT(z, t=tax)
            @test_throws ArgumentError VTTrac.setup(vtt, 2, 2; ixhw=3, iyhw=3, peak_inside_th=0.03)
        end

        @testset "VTT validates mask size" begin
            badmask = falses(nt-1, ny, nx)
            @test_throws ArgumentError VTTrac.VTT(z, t=tax, mask=badmask)
        end

        @testset "VTT requires at least 2 time steps" begin
            z1 = Array{Float32,3}(z[1:1,:,:])
            @test_throws ArgumentError VTTrac.VTT(z1)
        end
    end
end
