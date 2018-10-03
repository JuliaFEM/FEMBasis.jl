[![Build Status](https://travis-ci.org/JuliaFEM/FEMBasis.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/FEMBasis.jl)[![Coverage Status](https://coveralls.io/repos/github/JuliaFEM/FEMBasis.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaFEM/FEMBasis.jl?branch=master)[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliafem.github.io/FEMBasis.jl/stable)[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliafem.github.io/FEMBasis.jl/latest)[![Issues](https://img.shields.io/github/issues/JuliaFEM/FEMBasis.jl.svg)](https://github.com/JuliaFEM/FEMBasis.jl/issues)

`FEMBasis.jl` contains interpolation routines for standard finite element
function spaces.  Given ansatz and coordinates of domain, interpolation
functions are calculated  symbolically in a very general way to get efficient
code. As a concrete example, to generate basis functions for a standard 10-node
tetrahedron one can write

```julia
code = FEMBasis.create_basis(
    :Tet10,
    "10 node quadratic tetrahedral element",
    [
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
    ],
    :(1 + u + v + w + u*v + v*w + w*u + u^2 + v^2 + w^2),
   )
```

The resulting code is
```julia
    struct Tet10 <: FEMBasis.AbstractBasis{3}
    end
    Base.@pure function Base.size(::Type{Tet10})
            return (3, 10)
        end
    function Base.size(::Type{Tet10}, j::Int)
        j == 1 && return 3
        j == 2 && return 10
    end
    Base.@pure function Base.length(::Type{Tet10})
            return 10
        end
    function FEMBasis.get_reference_element_coordinates(::Type{Tet10})
        return Tensors.Tensor{1,3,Float64,3}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
    end
    function FEMBasis.eval_basis!(::Type{Tet10}, N::Vector{<:Number}, xi::Vec)
        assert length(N) == 10
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
    function FEMBasis.eval_dbasis!(::Type{Tet10}, dN::Vector{<:Vec{3}}, xi::Vec)
        @assert length(dN) == 10
        (u, v, w) = xi
        begin
            dN[1] = Vec(-3.0 + 4.0v + 4.0w + 2.0 * (2u), -3.0 + 4.0u + 4.0w + 2.0 * (2v), -3.0 + 4.0v + 4.0u + 2.0 * (2w))
            dN[2] = Vec(-1 + 2.0 * (2u), 0, 0)
            dN[3] = Vec(0, -1 + 2.0 * (2v), 0)
            dN[4] = Vec(0, 0, -1 + 2.0 * (2w))
            dN[5] = Vec(4.0 + -4.0v + -4.0w + -4.0 * (2u), -4.0u, -4.0u)
            dN[6] = Vec(4.0v, 4.0u, 0)
            dN[7] = Vec(-4.0v, 4.0 + -4.0u + -4.0w + -4.0 * (2v), -4.0v)
            dN[8] = Vec(-4.0w, -4.0w, 4.0 + -4.0v + -4.0u + -4.0 * (2w))
            dN[9] = Vec(4.0w, 0, 4.0u)
            dN[10] = Vec(0, 4.0w, 4.0v)
        end
        return dN
    end
end
```

Also more unusual elements can be defined. For example, pyramid element cannot be
descibed with ansatz, but it's still possible to implement by defining shape functions,
`Calculus.jl` is taking care of defining partial derivatives of function:
```julia
code = FEMBasis.create_basis(
    :Pyr5,
    "5 node linear pyramid element",
    [
     (-1.0, -1.0, -1.0), # N1
     ( 1.0, -1.0, -1.0), # N2
     ( 1.0,  1.0, -1.0), # N3
     (-1.0,  1.0, -1.0), # N4
     ( 0.0,  0.0,  1.0), # N5
    ],
    [
     :(1/8 * (1-u) * (1-v) * (1-w)),
     :(1/8 * (1+u) * (1-v) * (1-w)),
     :(1/8 * (1+u) * (1+v) * (1-w)),
     :(1/8 * (1-u) * (1+v) * (1-w)),
     :(1/2 * (1+w)),
    ],
   )
eval(code)
```

Basis function can have internal variables if needed, e.g. variable dof basis like
hierarchical basis functions or NURBS.

It's also possible to do some very common FEM calculations, like calculate Jacobian
or gradient of some variable with respect to some coordinates. For example, to 
calculate displacement gradient du/dX in unit square [0,1]^2, one could write:

```julia
using Tensors
B = Quad4()
X = Vec.([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
u = Vec.([(0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0)])
grad(B, u, X, Vec(0.0, 0.0))
```

Result is
```julia
2×2 Tensors.Tensor{2,2,Float64,4}:
 1.5  0.5
 1.0  2.0
```
