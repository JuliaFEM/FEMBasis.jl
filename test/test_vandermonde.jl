# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using FEMBasis: vandermonde_matrix
using Test

polynomial = :(1 + u + v + u*v)
coordinates = [(-1.0,-1.0), (1.0,-1.0), (1.0,1.0), (-1.0,1.0)]
V = vandermonde_matrix(polynomial, coordinates)
@debug "create Vandermonde matrix" polynomial coordinates V
V_expected = [
    1.0 -1.0 -1.0  1.0
    1.0  1.0 -1.0 -1.0
    1.0  1.0  1.0  1.0
    1.0 -1.0  1.0 -1.0]
@test isapprox(V, V_expected)

polynomial = :(1 + u + v)
coordinates = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
V = vandermonde_matrix(polynomial, coordinates)
V_expected = [
    1.0 0.0 0.0
    1.0 1.0 0.0
    1.0 0.0 1.0]
@test isapprox(V, V_expected)
