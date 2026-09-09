# `make_veech_right_triangle_and_basis`/`make_veech_right_triangle` live in
# QuantumBilliards.jl/src/utils/billiardutils.jl (not exported, so they're
# accessed here as `QuantumBilliards.<name>`) — do not redefine them here.

# solver: Vergini-Saraceno
# basis: corner adapted Fourier-Bessel
# billiard: Triangle
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate, compute_psi
@testset "Decomposition Method - Veech Triangle - Ground State" begin
    billiard, basis = QuantumBilliards.make_veech_right_triangle_and_basis(5)
    dim_scaling_factor = 2.0
    pts_scaling_factor = 5.0
    solver = VerginiSaracenoSolver(dim_scaling_factor, pts_scaling_factor)
    k0 = 6.1
    dk = 0.1
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 6.065082959967892
    t1_test = 0.0024524576991433386
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1912543992438375, 0.38049577139821616, 0.5657412040516646, 0.7450674834409242, 0.0, 0.38013047151484103, 0.7562509934152326, 1.1244116847521282, 1.4807817426748695, 0.0, 0.5642883950409858, 1.1226019179633464, 1.6690553994260682, 2.1979414922427374, 0.0, 0.741464683440748, 1.475037052816142, 2.1929436548777916, 2.887646656952914]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

@testset "Particular Solutions Method - Veech Triangle - Ground State" begin
    billiard, basis = QuantumBilliards.make_veech_right_triangle_and_basis(5)
    dim_scaling_factor = 2.0
    pts_scaling_factor = 5.0
    int_pts_scaling_factor = 2.0
    solver = ParticularSolutionsMethod(dim_scaling_factor, pts_scaling_factor, int_pts_scaling_factor)
    k0 = 6.1
    dk = 0.1
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 6.065082959967892
    t1_test = 3.117209745112226e-5
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0011315632476348654, 0.0022512163155965046, 0.0033472271992600958, 0.004408217104218537, 0.0, 0.0022490550090337384, 0.004474384942960184, 0.006652620125943378, 0.008761095413847607, 0.0, 0.003338631615866006, 0.006641912592537861, 0.009875023021577849, 0.013004195088027443, 0.0, 0.004386901041795259, 0.00872710680400177, 0.01297462525230263, 0.01708485841221074]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

 
# solver: Vergini-Saraceno
# basis: corner adapted Fourier-Bessel
# billiard: Triangle
# symmetry: None
# functions to test: solve_wavenumber, solve_spectrum, compute_eigenstate, compute_psi
@testset "Vergini Saraceno - Veech Triangle - Low Spectrum" begin
    billiard, basis = QuantumBilliards.make_veech_right_triangle_and_basis(5)
    dim_scaling_factor = 2.0
    pts_scaling_factor = 5.0
    solver = VerginiSaracenoSolver(dim_scaling_factor, pts_scaling_factor)
    k0 = 110.0
    dk = 0.1
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    ks, tens = solve_spectrum(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 110.0894534647971
    t1_test = 0.015990843403216694
    ks_test = [110.0894534647971]
    tens_test = [0.015990843403216694]
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, -5.143342048560577e-17, -0.07369796445157838, -0.5759107060556197, 0.485092827858497, -0.07082690430468566, 8.209052853850316e-16, 0.6633364271254106, 0.012220905610909936, 0.17194219594285254, -0.394705847179626, -5.337033860983294e-16, -0.5158790845472846, -0.13781489604910854, -0.26707812142865, 0.641073640804324, -1.6269090815529952e-15, -0.17088610978445745, 0.8049670084279076, -0.29475556993383234, -0.5235502340620063]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

# solver: Vergini-Saraceno
# basis: corner adapted Fourier-Bessel
# billiard: Triangle
# symmetry: None
# functions to test: solve_wavenumber, solve_spectrum, compute_eigenstate, compute_psi
@testset "Vergini Saraceno - Veech Triangle - High Spectrum" begin
    billiard, basis = QuantumBilliards.make_veech_right_triangle_and_basis(5)
    dim_scaling_factor = 2.0
    pts_scaling_factor = 5.0
    solver = VerginiSaracenoSolver(dim_scaling_factor, pts_scaling_factor)
    k0 = 2010.0
    dk = 0.05
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    ks, tens = solve_spectrum(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 2009.976791687898
    t1_test = 0.0010772639398423669
    ks_test = [2009.976791687898, 2010.0440083055544]
    tens_test = [0.0010772639398423669, 0.003873377109687749]
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 3.153419303678921e-16, -0.03218069856405793, -0.05992443758375497, 0.02218599885855573, 0.008485913632012934, -1.4988779039986679e-15, 0.029885022297499585, 0.07572240764904509, -0.04779010646702746, -0.018425717058992836, -5.274123942787825e-16, -0.016199467427262703, -0.08879585787572551, 0.04848742375577114, -0.026326463281711815, 4.696764947734056e-15, -0.03436979264720319, 0.03440559968103459, -0.07352007011994298, 0.0478196497413564]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

# solver: Vergini-Saraceno
# basis: Real Plane Waves
# billiard: Stadium
# symmetry: odd-odd
# functions to test: solve_wavenumber, solve_spectrum, compute_eigenstate, compute_psi
@testset "Vergini Saraceno - Stadium - Low Spectrum" begin
    billiard = StadiumBilliard(0.5)
    basis = RealPlaneWaves(12, sym_x = -1, sym_y = -1)
    dim_scaling_factor = 5.0
    pts_scaling_factor = 10.0
    solver = VerginiSaracenoSolver(dim_scaling_factor, pts_scaling_factor)
    k0 = 110.0
    dk = 0.05
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    ks, tens = solve_spectrum(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 109.98398490747392
    t1_test = 0.0005130410745052519
    ks_test = [109.98398490747392, 110.04189086784442]
    tens_test = [0.0005130410745052519, 0.0035083536721268477]
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04843409282971689, 0.09665875801932651, -0.05977820046546925, 0.014070514834491898, 0.0, -0.05562992330604402, 0.05239834305716773, -0.12039539147143066, 0.17401509881437177, 0.0, -0.008851265467139505, 0.10416057362570547, -0.12410238290507158, 0.076036687764999, 0.0, -0.0742024623925958, 0.044574433545681996, -0.0035890207209023917, -0.0070418588013708705]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

# solver: Vergini-Saraceno
# basis: Real Plane Waves
# billiard: Stadium
# symmetry: odd-odd
# functions to test: solve_wavenumber, solve_spectrum, compute_eigenstate, compute_psi
@testset "Vergini Saraceno - Stadium - High Spectrum" begin
    billiard = StadiumBilliard(0.5)
    basis = RealPlaneWaves(12, sym_x = -1, sym_y = -1)
    dim_scaling_factor = 5.0
    pts_scaling_factor = 10.0
    solver = VerginiSaracenoSolver(dim_scaling_factor, pts_scaling_factor)
    k0 = 1010.0
    dk = 0.01
    k, t1 = solve_wavenumber(solver, basis, billiard, k0, dk)
    ks, tens = solve_spectrum(solver, basis, billiard, k0, dk)
    state = compute_eigenstate(solver, basis, billiard, k)
    x_grid = collect(range(0.0, 0.1, length=5))
    y_grid = collect(range(0.0, 0.1, length=5))
    Psi = compute_psi(state, x_grid, y_grid; inside_only=true, memory_limit = 2.0e9, multithreaded = true)

    k_test = 1010.0004007793259
    t1_test = 3.212480087718333e-7
    ks_test = [1009.9936149553375, 1009.9961297133183, 1010.0004007793259, 1010.0015160257695]
    tens_test = [8.153810615325211e-5, 2.995835279701764e-5, 3.212480087718333e-7, 4.596661368602789e-6]
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02230935245832136, -0.04827226773775997, -0.012212416773368137, 0.08171119522093628, 0.0, 0.09992477276356589, 0.11625971043572127, -0.0662783727227309, 0.05570654869870634, 0.0, -0.11886122366282517, -0.03423830387400176, -0.14286513515755456, 0.07107632009292175, 0.0, 0.09245609513565874, -0.01263425673521082, -0.09868118878446691, 0.007432480267007546] 
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end

# solver: Double Layer Potential (Kress-corrected boundary integral method)
# basis: None (boundary-integral density, no basis expansion)
# billiard: Triangle (full, un-reduced boundary, all edges SpecularReflection)
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate
@testset "Double Layer Potential - Full Triangle - Ground State" begin
    billiard = QuantumBilliards.make_veech_right_triangle(5)
    pts_scaling_factor = 5.0
    solver = DoubleLayerPotentialSolver(pts_scaling_factor; grading=GlobalCornerGrading())
    k0 = 6.1
    dk = 0.1
    k, t1 = solve_wavenumber(solver, billiard, k0, dk)
    state = compute_eigenstate(solver, billiard, k)

    k_test = 6.065090955664264
    t1_test = 1.9925023646322987e-8
    ten_test = 1.992502365850128e-8
    dim_test = 200
    atol = 1e-3
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test isapprox(real(state.k), k_test; atol=atol)
    @test isapprox(state.ten, ten_test; atol=atol)
    @test state.dim == dim_test
    @test length(state.vec) == dim_test
end

# solver: Double Layer Potential (ungraded periodic boundary integral method)
# basis: None (boundary-integral density, no basis expansion)
# billiard: Circle (PolarBilliard, no true corners)
# symmetry: None
# functions to test: solve_wavenumber, k_sweep, compute_eigenstate
@testset "Double Layer Potential - Circle - Ground State" begin
    billiard = BilliardGeometry.PolarBilliard([0.0, 0.0])
    pts_scaling_factor = 5.0
    solver = DoubleLayerPotentialSolver(pts_scaling_factor; grading=SmoothPeriodicGrading())
    k0 = 2.4
    dk = 0.2
    k, t1 = solve_wavenumber(solver, billiard, k0, dk)
    ks = collect(range(2.35, 2.45, length=11))
    tens = k_sweep(solver, billiard, ks)
    state = compute_eigenstate(solver, billiard, k)

    k_test = 2.4048255557243405 # matches the first zero of J0, ≈ 2.404825557695772
    t1_test = 4.092042747328376e-9
    ks_test = [2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45]
    tens_test = [0.11387428129769038, 0.09309964793161138, 0.07232460472834643, 0.05155117343523519, 0.030781376045685045, 0.010017234588134472, 0.010739229084787815, 0.031485993506856644, 0.05222103780861643, 0.07294234192746765, 0.09364788681749107]
    ten_test = 4.0920427729930095e-9
    dim_test = 200
    atol = 1e-3
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test isapprox(state.ten, ten_test; atol=atol)
    @test state.dim == dim_test
    @test argmin(tens) == 6 # k=2.40 is the tension minimum in the swept window
end

# solver: Combined Field Integral Equation (Kress-corrected boundary integral method)
# basis: None (boundary-integral density, no basis expansion)
# billiard: Triangle (full, un-reduced boundary, all edges SpecularReflection)
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate
@testset "Combined Field Integral Equation - Full Triangle - Ground State" begin
    billiard = QuantumBilliards.make_veech_right_triangle(5)
    pts_scaling_factor = 5.0
    solver = CombinedFieldIntegralEquationSolver(pts_scaling_factor; grading=GlobalCornerGrading())
    k0 = 6.1
    dk = 0.1
    k, t1 = solve_wavenumber(solver, billiard, k0, dk)
    state = compute_eigenstate(solver, billiard, k)

    k_test = 6.065091021584176
    t1_test = 4.236448445773685e-8
    ten_test = 4.236448435637523e-8
    dim_test = 200
    atol = 1e-3
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test isapprox(real(state.k), k_test; atol=atol)
    @test isapprox(state.ten, ten_test; atol=atol)
    @test state.dim == dim_test
    @test length(state.vec) == dim_test
end

# solver: Combined Field Integral Equation (ungraded periodic boundary integral method)
# basis: None (boundary-integral density, no basis expansion)
# billiard: Circle (PolarBilliard, no true corners)
# symmetry: None
# functions to test: solve_wavenumber, k_sweep, compute_eigenstate
@testset "Combined Field Integral Equation - Circle - Ground State" begin
    billiard = BilliardGeometry.PolarBilliard([0.0, 0.0])
    pts_scaling_factor = 5.0
    solver = CombinedFieldIntegralEquationSolver(pts_scaling_factor; grading=SmoothPeriodicGrading())
    k0 = 2.4
    dk = 0.2
    k, t1 = solve_wavenumber(solver, billiard, k0, dk)
    ks = collect(range(2.35, 2.45, length=11))
    tens = k_sweep(solver, billiard, ks)
    state = compute_eigenstate(solver, billiard, k)

    k_test = 2.4048255524123827 # matches the first zero of J0, ≈ 2.404825557695772
    t1_test = 2.1430713503075886e-8
    ks_test = [2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45]
    tens_test = [0.2223184861379049, 0.1817923655089925, 0.14125066627230293, 0.10069757153072129, 0.060137263312625436, 0.01957392217808703, 0.020988273175142097, 0.06154514630575029, 0.10209252342310988, 0.14262623378153413, 0.18314211007456732]
    ten_test = 2.1430713549767833e-8
    dim_test = 200
    atol = 1e-3
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(ks, ks_test; atol=atol))
    @test all(isapprox.(tens, tens_test; atol=atol))
    @test isapprox(state.ten, ten_test; atol=atol)
    @test state.dim == dim_test
    @test argmin(tens) == 6 # k=2.40 is the tension minimum in the swept window
end

# solver: Double Layer Potential (Kress-corrected boundary integral method)
# basis: None (boundary-integral density, no basis expansion)
# billiard: Stadium (D2-symmetric quarter fundamental domain)
# symmetry: YAxisReflection (folds the complete physical boundary from
#           BilliardGeometry.full_boundary onto the fundamental domain)
# functions to test: evaluate_points, boundary_matrix_size, construct_matrices
#
# Regression test for the Step-4 migration-plan fix: evaluate_points used to
# discretize only the fundamental domain's own quarter boundary
# (get_boundary_curves) even when a symmetry was set, which is far too short
# a boundary for symmetry_index_orbits' exact index-permutation folding to be
# meaningful. It now discretizes the complete physical boundary
# (BilliardGeometry.full_boundary) whenever solver.symmetry !== nothing. This
# is verified by an exact algebraic identity rather than by locating a
# spectral resonance (which is sensitive to point-count/tolerance choices):
# the symmetry-reduced Fredholm matrix must equal the orbit-summed columns of
# the true (unreduced) Fredholm matrix assembled on the same full-boundary
# discretization.
@testset "Double Layer Potential - Stadium (YAxisReflection) - symmetry-reduced matrix consistency" begin
    billiard = StadiumBilliard(0.3)
    k = 5.5
    solver = DoubleLayerPotentialSolver(10.0; symmetry=BilliardGeometry.YAxisReflection())
    pts = evaluate_points(solver, billiard, k)

    # evaluate_points now samples the complete physical boundary, not just
    # the quarter fundamental domain.
    @test length(pts) > 4*length(BilliardGeometry.get_boundary_curves(billiard))

    A_reduced = construct_matrices(solver, pts, k)
    m = size(A_reduced, 1)
    @test m == boundary_matrix_size(solver, pts)

    graded = QuantumBilliards._is_nontrivial_dlp_grading(pts)
    G = QuantumBilliards.boundary_geom_cache(pts, graded)
    N = length(pts)
    Rmat = zeros(Float64, N, N)
    QuantumBilliards.kress_R!(Rmat)
    A_full = Matrix{ComplexF64}(undef, N, N)
    QuantumBilliards._dlp_fredholm_full!(A_full, pts, Rmat, G, k)

    orbits = BilliardGeometry.symmetry_index_orbits(Float64, pts.xy, solver.symmetry)
    fund = orbits.fundamental_indices
    expected = zeros(ComplexF64, m, m)
    for b in 1:m
        cols = findall(==(b), orbits.orbit_of)
        for a in 1:m
            expected[a,b] = sum(A_full[fund[a], j] for j in cols)
        end
    end
    @test maximum(abs.(expected .- A_reduced)) < 1e-10
end

