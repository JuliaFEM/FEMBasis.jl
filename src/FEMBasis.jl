# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

module FEMBasis

import Calculus
using LinearAlgebra
using Tensors
# Reexport Vec
export Vec

const Vecish{N, T} = Union{NTuple{N, T}, Vec{N, T}}

abstract type AbstractBasis{dim} end
# Forward methods on instances to types
Base.length(B::T) where T<:AbstractBasis = length(T)
Base.size(B::T) where T<:AbstractBasis = size(T)
eval_basis!(B::T, N, xi) where T<:AbstractBasis = eval_basis!(T, N, xi)
eval_dbasis!(B::T, dN, xi) where T<:AbstractBasis = eval_dbasis!(T, dN, xi)
eval_basis(B::AbstractBasis{dim}, xi) where {dim} = eval_basis!(B, zeros(length(B)), xi)
eval_dbasis(B::AbstractBasis{dim}, xi) where {dim} = eval_dbasis!(B, zeros(Vec{dim}, length(B)), xi)

include("subs.jl")
include("vandermonde.jl")

include("create_basis.jl")
export get_reference_element_coordinates, eval_basis!, eval_dbasis!

include("lagrange_segments.jl")
export Seg2, Seg3
include("lagrange_quadrangles.jl")
export Quad4, Quad8, Quad9
#=
include("lagrange_triangles.jl")
export Tri3, Tri6, Tri7
include("lagrange_tetrahedrons.jl")
export Tet4, Tet10
include("lagrange_hexahedrons.jl")
export Hex8, Hex20, Hex27
include("lagrange_wedges.jl")
export Wedge6, Wedge15
include("lagrange_pyramids.jl")
export Pyr5

include("nurbs.jl")
include("nurbs_segment.jl")
export NSeg
include("nurbs_surface.jl")
export NSurf
include("nurbs_solid.jl")
export NSolid
=#
include("math.jl")
export interpolate, interpolate!, jacobian
export grad, grad!
export BasisInfo

end
