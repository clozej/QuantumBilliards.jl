#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

"""
This struct represents a collection of parameters describing the Basis State in the default constructor scheme.

# Fields
- `k`: The wavenumber of the basis state.
- `k_basis`: The wavenumber of the basis state.
- `vec`: The coefficients of the basis state in the basis.
- `idx`: The index of the basis state in the basis.
- `dim`: The dimension of the basis.
- `eps`: The numerical precision of the underlying float type
- `basis`: The basis structure <: AbsState.
"""
struct BasisState{K,T,Ba} <: StationaryState 
    k::K
    k_basis::K
    vec::Vector{T}
    idx::Int64
    dim::Int64
    eps::T
    basis::Ba
end

"""
A constructor for a subtype of `AbsBasis` struct. The `dim` is iherited from the `basis`, the `eps` is also defined from the underlying float type of the `basis`.

# Arguments:
- `basis`: An instance of `AbsState` struct
- `k`: The wavenumber of the basis state.
- `i`: The index of the basis state in the basis.

# Returns:
- An instance of `AbsBasis` struct.
"""
function BasisState(basis, k, i)  
    dim = basis.dim
    typ = typeof(k)
    eps = set_precision(k)
    vec = zeros(typ,dim)
    vec[i] = one(typ)
    return BasisState(k,k, vec, i, dim, eps, basis)
end
