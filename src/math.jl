# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

"""
    interpolate(B, T, xi)

Given basis B, interpolate T at xi.

# Example
```jldoctest
B = Quad4()
X = Vec.([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
T = [1.0, 2.0, 3.0, 4.0]
interpolate(B, T, (0.0, 0.0))

# output

2.5
```
"""
function interpolate(B::AbstractBasis{dim}, T::Vector{<:Number}, xi::Vecish{dim}) where {dim}
    N = eval_basis(B, xi)
    return sum(b*t for (b, t) in zip(N, T))
end

"""
    jacobian(B, X, xi)

Given basis B, calculate jacobian at xi.

# Example
```jldoctest
B = Quad4()
X = Vec.([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
jacobian(B, X, Vec((0.0, 0.0)))

# output

2×2 Array{Float64,2}:
 0.5  0.0
 0.0  0.5

```
"""
jacobian(B::AbstractBasis{dim}, X, xi) where {dim} = jacobian(B, X, xi, eval_dbasis(B, xi))

function jacobian(B::AbstractBasis{dim}, X::Vector, xi, dB) where {dim}
    @assert length(X) == length(B)
    J = zero(Tensor{2, dim})
    for i = 1:length(X)
        J += otimes(dB[i], X[i]) # dB[i] ⊗ X[i]
    end
    return J
end



"""
    grad(B, X, xi)

Given basis B, calculate gradient dB/dX at xi.

# Example
```jldoctest
B = Quad4()
X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
grad(B, X, (0.0, 0.0))

# output

2×4 Array{Float64,2}:
 -0.5   0.5  0.5  -0.5
 -0.5  -0.5  0.5   0.5

```
"""
grad(B, X, xi) = _grad(B, X, xi, eval_dbasis(B, xi))

function _grad(B, X, xi, dB::Vector{<:Vec})
    J = jacobian(B, X, xi, dB)
    return [inv(J) ⋅ db for db in dB]
end

"""
    grad(B, T, X, xi)

Calculate gradient of `T` with respect to `X` in point `xi` using basis `B`.

# Example
```jldoctest
B = Quad4()
X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
u = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])
grad(B, u, X, (0.0, 0.0))

# output

2×2 Array{Float64,2}:
 1.5  0.5
 1.0  2.0

```
"""
function grad(B, T, X, xi)
    G = grad(B, X, xi)
    dTdX = sum(T[i] ⊗ G[i] for i=1:nbasis)
    return dTdX
end


"""
Data type for fast FEM.
"""
mutable struct BasisInfo{B<:AbstractBasis,dim, T}
    N::Vector{T}
    dN::Vector{T}
    grad::Vector{T}
    J::Tensor{2, dim, T}
    invJ::Tensor{2, dim, T}
    detJ::T
end

Base.length(B::BasisInfo{T}) where T<:AbstractBasis = length(T)
Base.size(B::BasisInfo{T})   where T<:AbstractBasis = size(T)

"""
Initialization of data type `BasisInfo`.

# Examples

```jldoctest

BasisInfo(Tri3)

# output

FEMBasis.BasisInfo{FEMBasis.Tri3,Float64}([0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0], 0.0)

```

"""
function BasisInfo(::Type{B}, T=Float64) where B <: AbstractBasis{dim} where dim
    nbasis = length(B)
    N = zeros(T, nbasis)
    dN = zeros(Vec{dim, T}, nbasis)
    grad = zeros(Vec{dim, T}, nbasis)
    J = zero(Tensor{2, dim, T})
    invJ = zero(Tensor{2, dim, T})
    detJ = zero(T)
    return BasisInfo{B,dim,T}(N, dN, grad, J, invJ, detJ)
end

"""
Evaluate basis, gradient and so on for some point `xi`.

# Examples

```jldoctest

b = BasisInfo(Quad4)
X = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))
xi = (0.0, 0.0)
eval_basis!(b, X, xi)

# output

FEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)

```
"""
function eval_basis!(bi::BasisInfo{B}, X, xi) where B
    # evaluate basis and derivatives
    eval_basis!(B, bi.N, xi)
    eval_dbasis!(B, bi.dN, xi)

    # calculate Jacobian
    bi.J = jacobian(B, X, xi, bi.dN)

    # calculate determinant of Jacobian + gradient operator

    # TODO, fixup curve + manifold
    @assert dims[1] == dims[2]
    bi.inv = inv(J)
    bi.grad = bi.invJ * bi.dN
    mul!(bi.grad, bi.invJ, bi.dN)
    bi.detJ = det(J)
    #=
    elseif dim1 == 1 # curve
        bi.detJ = norm(bi.J)
    elseif dim1 == 2 # manifold
        bi.detJ = norm(cross(bi.J[1,:], bi.J[2,:]))
    end
    =#

    return bi
end

"""
    grad!(bi, gradu, u)

Evalute gradient ∂u/∂X and store result to matrix `gradu`. It is assumed
that `eval_basis!` has been already run to `bi` so it already contains
all necessary matrices evaluated with some `X` and `xi`.

# Example

First setup and evaluate basis using `eval_basis!`:
```jldoctest ex1
B = BasisInfo(Quad4)
X = Vec.([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)])
xi = Vec(0.0, 0.0)
eval_basis!(B, X, xi)

# output

FEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)

```

Next, calculate gradient of `u`:
```jldoctest ex1
u = Vec.([(0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0)])
grad(B, u)

# output

2×2 Array{Float64,2}:
 1.5  0.5
 1.0  2.0
```
"""
function grad(bi::BasisInfo{B} where B <: AbstractBasis{dim}, u)  where dim
    nbasis = length(B)
    gradu = zero(Tensor{2, B})
    for k=1:nbasis
        gradu += bi.grad * u[k]
    end
    return gradu
end
