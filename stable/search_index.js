var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "(Image: Build Status)(Image: Coverage Status)(Image: )(Image: )(Image: Issues)FEMBasis.jl contains interpolation routines for standard finite element function spaces.  Given ansatz and coordinates of domain, interpolation functions are calculated  symbolically in a very general way to get efficient code. As a concrete example, to generate basis functions for a standard 10-node tetrahedron one can writecode = create_basis(\n    :Tet10,\n    \"10 node quadratic tetrahedral element\",\n    (\n     (0.0, 0.0, 0.0), # N1\n     (1.0, 0.0, 0.0), # N2\n     (0.0, 1.0, 0.0), # N3\n     (0.0, 0.0, 1.0), # N4\n     (0.5, 0.0, 0.0), # N5\n     (0.5, 0.5, 0.0), # N6\n     (0.0, 0.5, 0.0), # N7\n     (0.0, 0.0, 0.5), # N8\n     (0.5, 0.0, 0.5), # N9\n     (0.0, 0.5, 0.5), # N10\n    ),\n    :(1 + u + v + w + u*v + v*w + w*u + u^2 + v^2 + w^2),\n   )The resulting code isbegin\n    mutable struct Tet10\n    end\n    function Base.size(::Type{Tet10})\n        return (3, 10)\n    end\n    function Base.size(::Type{Tet10}, j::Int)\n        j == 1 && return 3\n        j == 2 && return 10\n    end\n    function Base.length(::Type{Tet10})\n        return 10\n    end\n    function FEMBasis.get_reference_element_coordinates(::Type{Tet10})\n        return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))\n    end\n    function FEMBasis.eval_basis!{T}(::Type{Tet10}, N::Matrix{T}, xi::Tuple{T, T, T})\n        (u, v, w) = xi\n        begin\n            N[1] = 1.0 + -3.0u + -3.0v + -3.0w + 4.0 * (u * v) + 4.0 * (v * w) + 4.0 * (w * u) + 2.0 * u ^ 2 + 2.0 * v ^ 2 + 2.0 * w ^ 2\n            N[2] = -u + 2.0 * u ^ 2\n            N[3] = -v + 2.0 * v ^ 2\n            N[4] = -w + 2.0 * w ^ 2\n            N[5] = 4.0u + -4.0 * (u * v) + -4.0 * (w * u) + -4.0 * u ^ 2\n            N[6] = +(4.0 * (u * v))\n            N[7] = 4.0v + -4.0 * (u * v) + -4.0 * (v * w) + -4.0 * v ^ 2\n            N[8] = 4.0w + -4.0 * (v * w) + -4.0 * (w * u) + -4.0 * w ^ 2\n            N[9] = +(4.0 * (w * u))\n            N[10] = +(4.0 * (v * w))\n        end\n        return N\n    end\n    function FEMBasis.eval_dbasis!{T}(::Type{Tet10}, dN::Matrix{T}, xi::Tuple{T, T, T})\n        (u, v, w) = xi\n        begin\n            dN[1, 1] = -3.0 + 4.0v + 4.0w + 2.0 * (2u)\n            dN[2, 1] = -3.0 + 4.0u + 4.0w + 2.0 * (2v)\n            dN[3, 1] = -3.0 + 4.0v + 4.0u + 2.0 * (2w)\n            dN[1, 2] = -1 + 2.0 * (2u)\n            dN[2, 2] = 0\n            dN[3, 2] = 0\n            dN[1, 3] = 0\n            dN[2, 3] = -1 + 2.0 * (2v)\n            dN[3, 3] = 0\n            dN[1, 4] = 0\n            dN[2, 4] = 0\n            dN[3, 4] = -1 + 2.0 * (2w)\n            dN[1, 5] = 4.0 + -4.0v + -4.0w + -4.0 * (2u)\n            dN[2, 5] = -4.0u\n            dN[3, 5] = -4.0u\n            dN[1, 6] = 4.0v\n            dN[2, 6] = 4.0u\n            dN[3, 6] = 0\n            dN[1, 7] = -4.0v\n            dN[2, 7] = 4.0 + -4.0u + -4.0w + -4.0 * (2v)\n            dN[3, 7] = -4.0v\n            dN[1, 8] = -4.0w\n            dN[2, 8] = -4.0w\n            dN[3, 8] = 4.0 + -4.0v + -4.0u + -4.0 * (2w)\n            dN[1, 9] = 4.0w\n            dN[2, 9] = 0\n            dN[3, 9] = 4.0u\n            dN[1, 10] = 0\n            dN[2, 10] = 4.0w\n            dN[3, 10] = 4.0v\n        end\n        return dN\n    end\nendAlso more unusual elements can be defined. For example, pyramid element cannot be descibed with ansatz, but it's still possible to implement by defining shape functions, Calculus.jl is taking care of defining partial derivatives of function:code = create_basis(\n    :Pyr5,\n    \"5 node linear pyramid element\",\n    (\n     (-1.0, -1.0, -1.0), # N1\n     ( 1.0, -1.0, -1.0), # N2\n     ( 1.0,  1.0, -1.0), # N3\n     (-1.0,  1.0, -1.0), # N4\n     ( 0.0,  0.0,  1.0), # N5\n    ),\n    (\n     :(1/8 * (1-u) * (1-v) * (1-w)),\n     :(1/8 * (1+u) * (1-v) * (1-w)),\n     :(1/8 * (1+u) * (1+v) * (1-w)),\n     :(1/8 * (1-u) * (1+v) * (1-w)),\n     :(1/2 * (1+w)),\n    ),\n   )\neval(code)Basis function can have internal variables if needed, e.g. variable dof basis like hierarchical basis functions or NURBS.It's also possible to do some very common FEM calculations, like calculate Jacobian or gradient of some variable with respect to some coordinates. For example, to  calculate displacement gradient du/dX in unit square [0,1]^2, one could write:B = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\nu = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])\ngrad(B, u, X, (0.0, 0.0))Result is2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0"
},

{
    "location": "api.html#FEMBasis.BasisInfo",
    "page": "API documentation",
    "title": "FEMBasis.BasisInfo",
    "category": "Type",
    "text": "Data type for fast FEM.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.BasisInfo-Union{Tuple{B}, Tuple{Type{B},Any}, Tuple{Type{B}}} where B<:FEMBasis.AbstractBasis",
    "page": "API documentation",
    "title": "FEMBasis.BasisInfo",
    "category": "Method",
    "text": "Initialization of data type BasisInfo.\n\nExamples\n\n\nBasisInfo(Tri3)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Tri3,Float64}([0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0], 0.0)\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.NSeg",
    "page": "API documentation",
    "title": "FEMBasis.NSeg",
    "category": "Type",
    "text": "NURBS segment. \n\n\n\n"
},

{
    "location": "api.html#FEMBasis.eval_basis!-Union{Tuple{B}, Tuple{FEMBasis.BasisInfo{B,T} where T,Any,Any}} where B",
    "page": "API documentation",
    "title": "FEMBasis.eval_basis!",
    "category": "Method",
    "text": "Evaluate basis, gradient and so on for some point xi.\n\nExamples\n\n\nb = BasisInfo(Quad4)\nX = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))\nxi = (0.0, 0.0)\neval_basis!(b, X, xi)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad!-Union{Tuple{B}, Tuple{FEMBasis.BasisInfo{B,T} where T,Any,Any}} where B",
    "page": "API documentation",
    "title": "FEMBasis.grad!",
    "category": "Method",
    "text": "grad!(bi, gradu, u)\n\nEvalute gradient ∂u/∂X and store result to matrix gradu. It is assumed that eval_basis! has been already run to bi so it already contains all necessary matrices evaluated with some X and xi.\n\nExample\n\nFirst setup and evaluate basis using eval_basis!:\n\nB = BasisInfo(Quad4)\nX = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))\nxi = (0.0, 0.0)\neval_basis!(B, X, xi)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)\n\n\nNext, calculate gradient of u:\n\nu = ((0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0))\ngradu = zeros(2, 2)\ngrad!(B, gradu, u)\n\n# output\n\n2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad-NTuple{4,Any}",
    "page": "API documentation",
    "title": "FEMBasis.grad",
    "category": "Method",
    "text": "grad(B, T, X, xi)\n\nCalculate gradient of T with respect to X in point xi using basis B.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\nu = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])\ngrad(B, u, X, (0.0, 0.0))\n\n# output\n\n2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad-Tuple{Any,Any,Any}",
    "page": "API documentation",
    "title": "FEMBasis.grad",
    "category": "Method",
    "text": "grad(B, X, xi)\n\nGiven basis B, calculate gradient dB/dX at xi.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\ngrad(B, X, (0.0, 0.0))\n\n# output\n\n2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Tuple{Any,Any,Any}",
    "page": "API documentation",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(B, T, xi)\n\nGiven basis B, interpolate T at xi.\n\nExample\n\nB = Quad4()\nX = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))\nT = (1.0, 2.0, 3.0, 4.0)\ninterpolate(B, T, (0.0, 0.0))\n\n# output\n\n2.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.jacobian-Tuple{Any,Any,Any}",
    "page": "API documentation",
    "title": "FEMBasis.jacobian",
    "category": "Method",
    "text": "jacobian(B, X, xi)\n\nGiven basis B, calculate jacobian at xi.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\njacobian(B, X, (0.0, 0.0))\n\n# output\n\n2×2 Array{Float64,2}:\n 0.5  0.0\n 0.0  0.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_basis_coefficients-Tuple{Expr,Tuple}",
    "page": "API documentation",
    "title": "FEMBasis.calculate_basis_coefficients",
    "category": "Method",
    "text": "function calculate_basis_coefficients(polynomial::String, coordinates::Vararg{Tuple})\n\nCalculate \"interpolate coefficient matrix\" for some polynomial p.\n\nExamples\n\nThat is, if we have polynomial p = 1 + u + v and coordinates (0,0), (1,0), (0,1), we find A*p such that first row is the first coordinate, second row is second coordinate and so on:\n\njulia> p = \"1 + u + w\"\njulia> X = ((0.0,0.0), (1.0,0.0), (0.0,1.0))\njulia> calculate_basis_coefficient(p, X)\n[1.0 0.0 0.0 # <-- p(0.0,0.0) = 1.0     = [1.0 0.0 0.0] * [1.0, u, v]\n 1.0 1.0 0.0 # <-- p(1.0,0.0) = 1.0 + u = [1.0 1.0 0.0] * [1.0, u, v]\n 1.0 0.0 1.0] # <- p(0.0,1.0) = 1.0 + v = [1.0 0.0 1.0] * [1.0, u, v]\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_interpolation_polynomial_derivatives-Tuple{Any,Any}",
    "page": "API documentation",
    "title": "FEMBasis.calculate_interpolation_polynomial_derivatives",
    "category": "Method",
    "text": "Calculate derivatives of basis functions with respect to parameters u, v, w.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_interpolation_polynomials-Tuple{Expr,Array{T,2} where T}",
    "page": "API documentation",
    "title": "FEMBasis.calculate_interpolation_polynomials",
    "category": "Method",
    "text": "\n\n"
},

{
    "location": "api.html#",
    "page": "API documentation",
    "title": "API documentation",
    "category": "page",
    "text": "DocTestSetup = quote\n    using FEMBasis\nendModules = [FEMBasis]"
},

]}
