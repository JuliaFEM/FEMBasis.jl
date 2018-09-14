# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

eval_basis(B::AbstractBasis, xi) = eval_basis!(B, zeros(length(B)), xi)
eval_dbasis(B::AbstractBasis, xi) = eval_dbasis!(B, zeros(Vec{dim(B)}, length(B)), xi)

"""
    interpolate(B, T, xi)

Given basis B, interpolate T at xi.

# Example
```jldoctest
B = Quad4()
X = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
T = (1.0, 2.0, 3.0, 4.0)
interpolate(B, T, (0.0, 0.0))

# output

2.5

```
"""
interpolate(B::AbstractBasis, T, xi) = interpolate(B, T, xi, eval_basis(B, xi))

function interpolate(B::AbstractBasis, T, xi, N)
    Ti = sum(b*t for (b, t) in zip(N, T))
    return Ti
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
jacobian(B::AbstractBasis, X, xi) = jacobian(B, X, xi, eval_dbasis(B, xi))

function jacobian(B::AbstractBasis, X::Vector, xi, dB)
    @assert length(X) == length(B)
    J = zero(Tensor{2, dim(B)})
    for i = 1:length(X)
        J += dB[i] ⊗ X[i]
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
grad(B, X, xi) = grad(B, X, xi, eval_dbasis(B, xi))

function grad(B, X, xi, dB)
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
mutable struct BasisInfo{B<:AbstractBasis,T}
    # TODO, Fixup
    N::Vector{T}
    dN::Vector{T}
    grad::Vector{T}
    J::Tensor{2, T}
    invJ::Tensor{2, T}
    detJ::T
end

function length(B::BasisInfo{T}) where T<:AbstractBasis
    return length(T)
end

function size(B::BasisInfo{T}) where T<:AbstractBasis
    return size(T)
end

"""
Initialization of data type `BasisInfo`.

# Examples

```jldoctest

BasisInfo(Tri3)

# output

FEMBasis.BasisInfo{FEMBasis.Tri3,Float64}([0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0], 0.0)

```

"""
function BasisInfo(::Type{B}, T=Float64) where B<:AbstractBasis
    dim, nbasis = size(B)
    N = zeros(T, 1, nbasis)
    dN = zeros(T, dim, nbasis)
    grad = zeros(T, dim, nbasis)
    J = zeros(T, dim, dim)
    invJ = zeros(T, dim, dim)
    detJ = zero(T)
    return BasisInfo{B,T}(N, dN, grad, J, invJ, detJ)
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
X = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))
xi = (0.0, 0.0)
eval_basis!(B, X, xi)

# output

FEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)

```

Next, calculate gradient of `u`:
```jldoctest ex1
u = ((0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0))
grad(B, gradu)

# output

2×2 Array{Float64,2}:
 1.5  0.5
 1.0  2.0
```
"""
function grad(bi::BasisInfo{B}, u) where B
    dim, nbasis = size(B)
    gradu = zero(Tensor{2, dim(B)})
    for k=1:nbasis
        gradu += bi.grad * u[k]
    end
    return gradu
end
