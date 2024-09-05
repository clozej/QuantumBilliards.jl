#include("../../abstracttypes.jl")
#include("../../utils/billiardutils.jl")
#include("decompositions.jl")
#include("../samplers.jl")
include("scalingmethod.jl")

"""
Solve for the wavenumber closest (in the interval dk) to the reference `k` using an accelerated solver (either ScalingMethodA or ScalingMethodB). This is a compendium of methods described in Alex Barnett's implementation of the Vergini Saraceno method. For further reading there is the link: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html

# Arguments
- `solver::AcceleratedSolver`: The solver used for the problem.
- `basis::AbsBasis`: The basis functions used in the problem.
- `billiard::AbsBilliard`: The billiard domain being analyzed.
- `k`: The reference wave number or frequency.
- `dk`: The tolerance for selecting eigenvalues close to `k`.

# Returns
- `ks[idx]`: The wavenumber closest to the reference `k`.
- `ts[idx]`: The corresponding tension value.

# Logic
1. **Dimension Calculation**: The problem's dimension is determined based on the billiard's length and the reference `k`, scaled by the solver's dimension scaling factor. The minimum dimension is a failsafe for smaller k values where the matrix F & Fk would be too small. In almost all cases we have the usual b (scaling factor) determined size
2. **Basis Resizing**: The basis functions are resized to match the calculated dimension from step 1.
3. **Point Evaluation**: Boundary points are evaluated for the billiard at the given wavenumber `k`. This is internally done with the the solver that under the hood to contains the sampler for the boundary points.
4. **Solve Eigenproblem**: The generalized eigenproblem is solved to find potential wavenumbers and tensions.
5. **Select Closest Wavenumber**: The wavenumber closest to `k` is selected and returned along with its corresponding tension (in the interval dk).
"""
function solve_wavenumber(solver::AcceleratedSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk)
        dim = max(solver.min_dim,round(Int, billiard.length*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, ts = solve(solver,new_basis,pts,k,dk)
    idx = findmin(abs.(ks.-k))[2]
    return ks[idx], ts[idx]
end

"""
Solve for the spectrum (in the interval dk) of the reference `k` using an accelerated solver (either ScalingMethodA or ScalingMethodB). This is a compendium of methods described in Alex Barnett's implementation of the Vergini Saraceno method. For further reading there is the link: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html

# Arguments
- `solver::AcceleratedSolver`: The solver used for the problem.
- `basis::AbsBasis`: The basis functions used in the problem.
- `billiard::AbsBilliard`: The billiard domain being analyzed.
- `k`: The reference wave number or frequency.
- `dk`: The tolerance for selecting eigenvalues close to `k`.

# Returns
- `ks`: A vector of wavenumbers within the specified tolerance of `k`.
- `ts`: A vector of corresponding tension values.

# Logic
1. **Dimension Calculation**: The problem's dimension is determined based on the billiard's length and the reference `k`, scaled by the solver's dimension scaling factor. The minimum dimension is a failsafe for smaller k values where the matrix F & Fk would be too small. In almost all cases we have the usual b (scaling factor) determined size
2. **Basis Resizing**: The basis functions are resized to match the calculated dimension from step 1.
3. **Point Evaluation**: Boundary points are evaluated for the billiard at the given wavenumber `k`. This is internally done with the the solver that under the hood to contains the sampler for the boundary points.
4. **Solve Eigenproblem**: The generalized eigenproblem is solved to find potential wavenumbers and tensions.
5. **Select Closest Wavenumber**: Return all the wavenumbers and corresponding tensions within the dk range.
"""
function solve_spectrum(solver::AcceleratedSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk)
    dim = max(solver.min_dim,round(Int, billiard.length*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, ts = solve(solver,new_basis,pts,k,dk)
    return ks, ts
end