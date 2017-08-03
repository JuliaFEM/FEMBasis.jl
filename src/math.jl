# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

import Base: length, size

function length{T<:AbstractBasis}(B::T)
    return length(T)
end

function size{T<:AbstractBasis}(B::T)
    return size(T)
end

function eval_basis!{T<:AbstractBasis}(B::T, N, xi)
    return eval_basis!(T, N, xi)
end

function eval_dbasis!{T<:AbstractBasis}(B::T, dN, xi)
    return eval_dbasis!(T, dN, xi)
end

function eval_basis(B, xi)
    N = zeros(1, length(B))
    eval_basis!(B, N, xi)
    return N
end

function eval_dbasis(B, xi)
    dN = zeros(size(B)...)
    eval_dbasis!(B, dN, xi)
    return dN
end

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
function interpolate(B, T, xi)
    Ti = sum(b*t for (b, t) in zip(eval_basis(B, xi), T))
    return Ti
end

"""
    jacobian(B, X, xi)

Given basis B, calculate jacobian at xi.

# Example
```jldoctest
B = Quad4()
X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
jacobian(B, X, (0.0, 0.0))

# output

2×2 Array{Float64,2}:
 0.5  0.0
 0.0  0.5

```
"""
function jacobian(B, X, xi)
    dB = eval_dbasis(B, xi)
    nbasis = length(B)
    J = sum(kron(dB[:,i], X[i]') for i=1:nbasis)
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
function grad(B, X, xi)
    J = jacobian(B, X, xi)
    dB = eval_dbasis(B, xi)
    G = inv(J) * dB
    return G
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
    nbasis = length(B)
    dTdX = sum(kron(G[:,i], T[i]') for i=1:nbasis)
    return dTdX'
end
