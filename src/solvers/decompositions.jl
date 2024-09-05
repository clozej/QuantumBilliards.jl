using LinearAlgebra
#=
function GESVDVALS(A, B; eps=1e-14)
    M = reduce(vcat, [A, B]) #concatenate columns
    Q, R =  qr(M)
    _ , sv_r, Vt_r = svd(R)
    mask = sv_r .> eps
    #println(mask)
    V1t = Vt_r[ : , mask]

    return svdvals((A * V1t), (B * V1t))
end
=#

"""
Compute the generalized eigenvalues and eigenvectors for the matrix pair (A, B). This is actually the logic found in Alex Barnett's thesis which describes how we get rid of the nullspace of the A (F) matrix. The idx variable gives us all the indexes of the eigenvalues that are bigger than some threshold ϵ and also get their eigenvectors. These eigenvectors are then stored in S. For this reduced subspace we then solve C^T * B * C (which is a symmetric matrix so we make sure with Symmetric call) (we solve for the reduced subspace). Finally we solve the reduced eigenproblem for this C^T * B * C

    normalisation: https://users.flatironinstitute.org/~ahb/thesis_html/node60.html : Truncating the singular generalized eigenproblem
    solver logic: https://users.flatironinstitute.org/~ahb/thesis_html/node58.html : Representation by a Helmholtz basis


# Arguments
- `A`: The first matrix.
- `B`: The second matrix.
- `eps`: Numerical tolerance to filter small eigenvalues of `A`. Defaults to `1e-15`.

# Returns
- `mu`: The generalized eigenvalues.
- `Z`: The matrix of eigenvectors.
- `C`: The transformation matrix applied to the eigenvectors.
"""
function generalized_eigen(A,B;eps=1e-15)
    d, S = eigen(A)
    idx = (d/maximum(d)) .> eps
    q = 1.0 ./ sqrt.(d[idx])
    #println(length(q))
    C = q' .* S[:,idx] 
    D = B * C
    #println(size(D))
    E = Symmetric(C' * D)
    #println(size(E))
    mu, Z = eigen(E) #check eigenvectors
    return mu, Z, C
end

"""
Compute the generalized eigenvalues and eigenvectors for the matrix pair (A, B). This is actually the logic found in Alex Barnett's thesis which describes how we get rid of the nullspace of the A (F) matrix. The idx variable gives us all the indexes of the eigenvalues that are bigger than some threshold ϵ and also get their eigenvectors. These eigenvectors are then stored in S. For this reduced subspace we then solve C^T * B * C (which is a symmetric matrix so we make sure with Symmetric call) (we solve for the reduced subspace). Finally we solve the reduced eigenproblem for this C^T * B * C

    normalisation: https://users.flatironinstitute.org/~ahb/thesis_html/node60.html : Truncating the singular generalized eigenproblem
    solver logic: https://users.flatironinstitute.org/~ahb/thesis_html/node58.html : Representation by a Helmholtz basis

# Arguments
- `A`: The first matrix.
- `B`: The second matrix.
- `eps`: Numerical tolerance to filter small eigenvalues of `A`. Defaults to `1e-15`.

# Returns
- `mu`: The generalized eigenvalues.
"""
function generalized_eigvals(A,B;eps=1e-15)
    d, S = eigen(A)
    idx = (d/maximum(d)) .> eps
    q = 1.0 ./ sqrt.(d[idx])
    #println(length(q))
    C = q' .* S[:,idx] 
    #println(size(C))
    D = B * C
    #println(size(D))
    E = Symmetric(C' * D)
    #println(size(E))
    mu = eigvals(E)
    return mu
end

"""
Compute the direct sum of two matrices.

# Arguments
- `A`: The first matrix.
- `B`: The second matrix.

# Returns
- A block diagonal matrix with `A` in the top-left block, `B` in the bottom-right block, and zeros in the off-diagonal blocks.

# Logic
The direct sum places matrix `A` in the top-left corner and matrix `B` in the bottom-right corner, with the remaining positions filled with zeros. The resulting matrix effectively combines `A` and `B` into a single larger matrix while keeping them in separate subspaces.
"""
directsum(A,B) = [A zeros(size(A,1), size(B,2)); zeros(size(B,1), size(A,2)) B]

"""
Adjust the scaling factors and samplers to match the number of boundary curves in a billiard problem.

# Arguments
- `solver::AbsSolver`: The solver instance containing the initial scaling factors and samplers.
- `billiard::AbsBilliard`: The billiard problem containing boundary curves.

# Returns
- `bs`: Adjusted vector of scaling factors, one for each boundary curve.
- `samplers`: Adjusted vector of samplers, one for each boundary curve.

# Logic
This function ensures that the number of scaling factors (`bs`) and samplers matches the number of curves in the billiard's fundamental boundary:

1. **Count Curves**: It first counts the number of curves in the billiard's boundary that are subtypes of `AbsRealCurve`.
2. **Adjust Scaling Factors (`bs`)**: If the number of curves exceeds the length of the `bs` vector, it extends `bs` by repeating the smallest scaling factor until it matches the number of curves.
3. **Adjust Samplers**: Similarly, if the number of curves exceeds the number of samplers, it extends the `samplers` vector by repeating the first sampler in the list until it matches the number of curves.

This adjustment ensures that every curve has a corresponding scaling factor and sampler, preventing mismatches during computations.
"""
function adjust_scaling_and_samplers(solver::AbsSolver, billiard::AbsBilliard)
    bs = solver.pts_scaling_factor
    samplers = solver.sampler
    default = samplers[1]
    n_curves = 0
    for crv in billiard.fundamental_boundary
        if typeof(crv) <: AbsRealCurve
            n_curves += 1
        end
    end

    b_min = minimum(bs)
    while length(bs)<n_curves
        push!(bs, b_min)
    end
    
    while length(samplers)<n_curves
        push!(samplers, default)
    end
    return bs, samplers
end
