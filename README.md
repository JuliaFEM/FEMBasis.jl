# FEMBasis.jl

[![Build Status](https://travis-ci.org/JuliaFEM/FEMBasis.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/FEMBasis.jl)[![Coverage Status](https://coveralls.io/repos/github/JuliaFEM/FEMBasis.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaFEM/FEMBasis.jl?branch=master)[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliafem.github.io/FEMBasis.jl/stable)[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliafem.github.io/FEMBasis.jl/latest)[![Issues](https://img.shields.io/github/issues/JuliaFEM/FEMBasis.jl.svg)](https://github.com/JuliaFEM/FEMBasis.jl/issues)

Package contains interpolation routines for standard finite element function spaces. 
Given ansatz and coordinates of domain, interpolation functions are calculated 
symbolically in a very general way to get efficient code. As a concrete example, 
to evaluate basis functions for standard tetrahedron we write

```julia
code = create_basis(
    :Tet10,
    "10 node quadratic tetrahedral element",
    (
     (0.0, 0.0, 0.0), # N1
     (1.0, 0.0, 0.0), # N2
     (0.0, 1.0, 0.0), # N3
     (0.0, 0.0, 1.0), # N4
     (0.5, 0.0, 0.0), # N5
     (0.5, 0.5, 0.0), # N6
     (0.0, 0.5, 0.0), # N7
     (0.0, 0.0, 0.5), # N8
     (0.5, 0.0, 0.5), # N9
     (0.0, 0.5, 0.5), # N10
    ),
    "1 + u + v + w + u*v + v*w + w*u + u^2 + v^2 + w^2",
   )
```

The resulting code is
```julia
begin
    mutable struct Tet10
    end
    function Base.size(::Type{Tet10})
        return (3, 10)
    end
    function Base.size(::Type{Tet10}, j::Int)
        j == 1 && return 3
        j == 2 && return 10
    end
    function Base.length(::Type{Tet10})
        return 10
    end
    function FEMBasis.get_reference_element_coordinates(::Type{Tet10})
        return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    end
    function FEMBasis.eval_basis!{T}(::Type{Tet10}, N::Matrix{T}, xi::Tuple{T, T, T})
        (u, v, w) = xi
        begin
            N[1] = 1.0 + -3.0u + -3.0v + -3.0w + 4.0 * (u * v) + 4.0 * (v * w) + 4.0 * (w * u) + 2.0 * u ^ 2 + 2.0 * v ^ 2 + 2.0 * w ^ 2
            N[2] = -u + 2.0 * u ^ 2
            N[3] = -v + 2.0 * v ^ 2
            N[4] = -w + 2.0 * w ^ 2
            N[5] = 4.0u + -4.0 * (u * v) + -4.0 * (w * u) + -4.0 * u ^ 2
            N[6] = +(4.0 * (u * v))
            N[7] = 4.0v + -4.0 * (u * v) + -4.0 * (v * w) + -4.0 * v ^ 2
            N[8] = 4.0w + -4.0 * (v * w) + -4.0 * (w * u) + -4.0 * w ^ 2
            N[9] = +(4.0 * (w * u))
            N[10] = +(4.0 * (v * w))
        end
        return N
    end
    function FEMBasis.eval_dbasis!{T}(::Type{Tet10}, dN::Matrix{T}, xi::Tuple{T, T, T})
        (u, v, w) = xi
        begin
            dN[1, 1] = -3.0 + 4.0v + 4.0w + 2.0 * (2u)
            dN[2, 1] = -3.0 + 4.0u + 4.0w + 2.0 * (2v)
            dN[3, 1] = -3.0 + 4.0v + 4.0u + 2.0 * (2w)
            dN[1, 2] = -1 + 2.0 * (2u)
            dN[2, 2] = 0
            dN[3, 2] = 0
            dN[1, 3] = 0
            dN[2, 3] = -1 + 2.0 * (2v)
            dN[3, 3] = 0
            dN[1, 4] = 0
            dN[2, 4] = 0
            dN[3, 4] = -1 + 2.0 * (2w)
            dN[1, 5] = 4.0 + -4.0v + -4.0w + -4.0 * (2u)
            dN[2, 5] = -4.0u
            dN[3, 5] = -4.0u
            dN[1, 6] = 4.0v
            dN[2, 6] = 4.0u
            dN[3, 6] = 0
            dN[1, 7] = -4.0v
            dN[2, 7] = 4.0 + -4.0u + -4.0w + -4.0 * (2v)
            dN[3, 7] = -4.0v
            dN[1, 8] = -4.0w
            dN[2, 8] = -4.0w
            dN[3, 8] = 4.0 + -4.0v + -4.0u + -4.0 * (2w)
            dN[1, 9] = 4.0w
            dN[2, 9] = 0
            dN[3, 9] = 4.0u
            dN[1, 10] = 0
            dN[2, 10] = 4.0w
            dN[3, 10] = 4.0v
        end
        return dN
    end
end
```

Also more unusual elements can be defined. For example, pyramid element cannot be
descibed with ansatz, but it's still possible to implement by defining shape functions,
`Calculus.jl` is taking care of defining partial derivatives of function:
```julia
code = create_basis(
    :Pyr5,
    "5 node linear pyramid element",
    (
     (-1.0, -1.0, -1.0), # N1
     ( 1.0, -1.0, -1.0), # N2
     ( 1.0,  1.0, -1.0), # N3
     (-1.0,  1.0, -1.0), # N4
     ( 0.0,  0.0,  1.0), # N5
    ),
    (
     "1/8 * (1-u) * (1-v) * (1-w)",
     "1/8 * (1+u) * (1-v) * (1-w)",
     "1/8 * (1+u) * (1+v) * (1-w)",
     "1/8 * (1-u) * (1+v) * (1-w)",
     "1/2 * (1+w)",
    ),
   )
eval(code)
```

Basis function can have internal variables if needed, e.g. variable dof basis like
hierarchical basis functions or NURBS.

It's also possible to do some very common FEM calculations, like calculate Jacobian
or gradient of some variable with respect to some coordinates. For example, to 
calculate displacement gradient du/dX in unit square [0,1]^2, one could write:

```julia
B = Quad4()
X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
u = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])
grad(B, u, X, (0.0, 0.0))
```

Result is
```julia
2×2 Array{Float64,2}:
 1.5  0.5
 1.0  2.0
```
