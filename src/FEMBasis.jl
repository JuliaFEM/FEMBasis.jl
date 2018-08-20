# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

module FEMBasis

using LinearAlgebra

abstract type AbstractBasis end

include("subs.jl")
include("vandermonde.jl")

include("create_basis.jl")
export get_reference_element_coordinates, eval_basis!, eval_dbasis!

include("lagrange_segments.jl")
export Seg2, Seg3
include("lagrange_quadrangles.jl")
export Quad4, Quad8, Quad9
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

include("math.jl")
export interpolate, interpolate!, jacobian
export grad, grad!
export BasisInfo

end
