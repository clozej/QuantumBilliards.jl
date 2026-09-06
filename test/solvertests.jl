
function make_veech_right_triangle_and_basis(n; edge_i=1)
    if n < 4
        println("Order must be 4 or above.")
        print("Setting n to 4.")
        n = 4
    end
    chi = (n-2)/2
    gamma = pi/2.0
    return make_triangle_and_basis(gamma, chi; edge_i)
end


# solver: Vergini-Saraceno
# basis: corner adapted Fourier-Bessel
# billiard: Triangle
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate, compute_psi
@testset "Decomposition Method - Veech Triangle - Ground State" begin
    billiard, basis = make_veech_right_triangle_and_basis(5)
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
#=
@testset "Particular Solutions Method - Veech Triangle - Ground State" begin
    billiard, basis = make_veech_right_triangle_and_basis(5)
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
    t1_test = 0.0024524576991433386
    psi_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1912543992438375, 0.38049577139821616, 0.5657412040516646, 0.7450674834409242, 0.0, 0.38013047151484103, 0.7562509934152326, 1.1244116847521282, 1.4807817426748695, 0.0, 0.5642883950409858, 1.1226019179633464, 1.6690553994260682, 2.1979414922427374, 0.0, 0.741464683440748, 1.475037052816142, 2.1929436548777916, 2.887646656952914]
    atol = 1e-3 
    @test isapprox(k, k_test; atol=atol)
    @test isapprox(t1, t1_test; atol=atol)
    @test all(isapprox.(Psi[:], psi_test; atol=atol))
end
=#
 
# solver: Vergini-Saraceno
# basis: corner adapted Fourier-Bessel
# billiard: Triangle
# symmetry: None
# functions to test: solve_wavenumber, solve_spectrum, compute_eigenstate, compute_psi
@testset "Vergini Saraceno - Veech Triangle - Low Spectrum" begin
    billiard, basis = make_veech_right_triangle_and_basis(5)
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
    billiard, basis = make_veech_right_triangle_and_basis(5)
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

