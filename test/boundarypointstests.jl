using StaticArrays
using LinearAlgebra

@testset "BoundaryPoints - extended fields and validation" begin
    T = Float64
    xy = [SVector{2,T}(cos(θ),sin(θ)) for θ in range(0.0,2pi;length=5)]
    # keyword constructor: unspecified fields default to empty/neutral
    bp = BoundaryPoints(xy)
    @test length(bp) == 5
    @test isempty(bp.normal) && isempty(bp.w) && isempty(bp.tangent)
    @test bp.compid == 1
    @test bp.is_periodic == true
    @test bp.shift_x == 0.0 && bp.shift_y == 0.0
    @test bp.xL == SVector(0.0,0.0)

    # length validation still enforced for both legacy and new fields
    @test_throws ErrorException BoundaryPoints(xy; normal=[SVector(1.0,0.0)])
    @test_throws ErrorException BoundaryPoints(xy; tangent=[SVector(1.0,0.0)])
    @test_throws ErrorException BoundaryPoints(xy; w=[1.0,2.0])

    # legacy fields (kappa, rdotn, w_vs, w_dm) remain unaffected
    bp2 = BoundaryPoints(xy; w_vs=fill(0.1,5), w_dm=fill(0.2,5))
    @test bp2.w_vs == fill(0.1,5)
    @test bp2.w_dm == fill(0.2,5)
end

@testset "BoundaryPoints - parametrized constructor (normal from tangent)" begin
    T = Float64
    N = 32
    R = 2.0
    ts = [2*pi*(k-0.5)/N for k in 1:N]
    tphys = ts ./ (2*pi)
    circ = BilliardGeometry.CircleSegment(R, 2*pi, 0.0, (0.0,0.0))
    xy = BilliardGeometry.curve(circ, tphys)
    tangent = BilliardGeometry.tangent(circ, tphys) ./ (2*pi)
    tangent_2 = BilliardGeometry.tangent_2(circ, tphys) ./ (2*pi)^2
    s = BilliardGeometry.arc_length(circ, tphys)
    h = 2*pi/N
    ds = [hypot(v[1],v[2])*h for v in tangent]
    ws = fill(h,N)
    ws_der = ones(N)
    z = SVector(0.0,0.0)
    bp = BoundaryPoints(xy, tangent, tangent_2, ts, tphys, ws, ws_der, s, ds, 1, true, z, z, z, z)

    @test length(bp) == N
    @test length(bp.normal) == N
    # hand-computed normal n = (t_y,-t_x)/|t| matches the constructor's computation
    for i in 1:N
        tx,ty = tangent[i]
        n_hand = SVector(ty/hypot(tx,ty), -tx/hypot(tx,ty))
        @test isapprox(bp.normal[i], n_hand; atol=1e-12)
    end
    # for a circle centered at the origin the outward normal is parallel to xy/|xy|
    for i in 1:N
        @test isapprox(bp.normal[i], xy[i]./norm(xy[i]); atol=1e-8)
    end
end

@testset "BoundaryPoints - multi-component helpers" begin
    T = Float64
    xy1 = [SVector{2,T}(Float64(i),0.0) for i in 1:3]
    xy2 = [SVector{2,T}(Float64(i),1.0) for i in 1:5]
    bp1 = BoundaryPoints(xy1; s=[0.0,1.0,2.0], ds=fill(1.0,3))
    bp2 = BoundaryPoints(xy2; s=[0.0,1.0,2.0,3.0,4.0], ds=fill(1.0,5))
    comps = [bp1,bp2]

    @test boundary_matrix_size(comps) == 8
    @test boundary_matrix_size(bp1) == 3
    @test component_offsets(comps) == [1,4,9]
    @test component_offsets(bp1) == [1,4]

    s_flat = boundary_s(comps)
    @test length(s_flat) == 8
    # second component's arc length is shifted by the first component's total length (3.0)
    @test isapprox(s_flat[4], 0.0+3.0; atol=1e-12)
    @test isapprox(s_flat[end], 4.0+3.0; atol=1e-12)
end

@testset "boundary_geom_cache / BoundaryPanelArrays" begin
    T = Float64
    N = 48
    R = 1.5
    ts = [2*pi*(k-0.5)/N for k in 1:N]
    tphys = ts ./ (2*pi)
    circ = BilliardGeometry.CircleSegment(R, 2*pi, 0.0, (0.0,0.0))
    xy = BilliardGeometry.curve(circ, tphys)
    tangent = BilliardGeometry.tangent(circ, tphys) ./ (2*pi)
    tangent_2 = BilliardGeometry.tangent_2(circ, tphys) ./ (2*pi)^2
    s = BilliardGeometry.arc_length(circ, tphys)
    h = 2*pi/N
    ds = [hypot(v[1],v[2])*h for v in tangent]
    ws = fill(h,N)
    ws_der = ones(N)
    z = SVector(0.0,0.0)
    bp = BoundaryPoints(xy, tangent, tangent_2, ts, tphys, ws, ws_der, s, ds, 1, true, z, z, z, z)

    panels = QuantumBilliards._boundary_panel_arrays_cache(bp)
    @test panels.X == getindex.(xy,1)
    @test panels.Y == getindex.(xy,2)
    @test all(isapprox.(panels.speed, R; atol=1e-8)) # |gamma'(t)| = R*2pi / (2pi) = R for this parametrization

    cache = boundary_geom_cache(bp)
    @test all(isfinite, cache.R)
    @test all(isfinite, cache.invR)
    @test all(cache.R[i,i] == one(T) for i in 1:N) # diagonal set to 1 before inversion, per convention
    @test all(cache.invR[i,i] == zero(T) for i in 1:N)
    @test isapprox(cache.R, cache.R'; atol=1e-10) # pairwise distance matrix is symmetric
    @test isempty(cache.original_ts) # corner_kress defaults to false

    cache_kress = boundary_geom_cache(bp, true)
    @test cache_kress.original_ts == ts

    # `cache.kappa` is `-(gamma'xgamma'')/|gamma'|^2 / (2pi)` (a dimensionless
    # quantity used in the diagonal Kress-kernel correction, not the raw
    # geometric curvature): for a circle parametrized by ts=angle this is the
    # constant `-1/(2pi)` independent of R, since curvature (1/R) is exactly
    # cancelled by the parametrization speed (R).
    @test all(isapprox.(cache.kappa, -1/(2*pi); atol=1e-6))

    nx,ny,speed = component_normals(bp)
    @test all(isapprox.(hypot.(nx,ny), 1.0; atol=1e-10)) # normals are unit vectors
    @test all(isapprox.(speed, R; atol=1e-8))
end

@testset "flatten_boundary_components / flatten_boundary_ds" begin
    T = Float64
    xy1 = [SVector{2,T}(0.0,0.0), SVector{2,T}(1.0,0.0)]
    xy2 = [SVector{2,T}(0.0,1.0), SVector{2,T}(1.0,1.0), SVector{2,T}(2.0,1.0)]
    n1 = [SVector{2,T}(0.0,-1.0), SVector{2,T}(0.0,-1.0)]
    n2 = [SVector{2,T}(0.0,1.0), SVector{2,T}(0.0,1.0), SVector{2,T}(0.0,1.0)]
    bp1 = BoundaryPoints(xy1; normal=n1, ds=fill(0.5,2))
    bp2 = BoundaryPoints(xy2; normal=n2, ds=fill(0.25,3))
    comps = [bp1,bp2]

    flat = flatten_boundary_components(comps)
    @test length(flat.x) == 5
    @test flat.x == [0.0,1.0,0.0,1.0,2.0]
    @test flat.y == [0.0,0.0,1.0,1.0,1.0]
    @test flat.ny == [-1.0,-1.0,1.0,1.0,1.0]
    @test flat.offs == [1,3,6]

    ds_flat = flatten_boundary_ds(comps)
    @test ds_flat == [0.5,0.5,0.25,0.25,0.25]
end

@testset "points_in_billiard" begin
    billiard, _ = make_veech_right_triangle_and_basis(4)
    interior = QuantumBilliards.random_interior_points(billiard, 5)
    @test all(points_in_billiard(interior, billiard))
    far_pt = SVector(1e3,1e3)
    @test points_in_billiard([far_pt], billiard) == [false]
end

@testset "estimate_rmin_rmax" begin
    xy = SVector{2,Float64}[(0,0),(1,0),(1,1),(0,1)]
    pts = BoundaryPoints(xy)
    rmin, rmax = estimate_rmin_rmax(pts, nothing)
    @test isapprox(rmin, 1.0; atol=1e-10)
    @test isapprox(rmax, sqrt(2); atol=1e-10)

    # symmetry-reduced variant, evaluated on the complete physical boundary
    # discretization (BilliardGeometry.full_boundary): sanity-checked against
    # the direct pairwise version on the same points.
    billiard = StadiumBilliard(0.3)
    solver = DoubleLayerPotentialSolver(10.0; symmetry=BilliardGeometry.YAxisReflection())
    bpts = evaluate_points(solver, billiard, 5.5)
    rmin_full, rmax_full = estimate_rmin_rmax(bpts, nothing)
    rmin_sym, rmax_sym = estimate_rmin_rmax(bpts, solver.symmetry)
    @test rmin_sym > 0 && isfinite(rmin_sym)
    @test rmax_sym > 0 && isfinite(rmax_sym)
    @test rmax_sym <= rmax_full + 1e-8
end
