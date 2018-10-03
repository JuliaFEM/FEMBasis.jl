# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using FEMBasis: calculate_interpolation_polynomials,
                calculate_interpolation_polynomial_derivatives

p = :(1 + u + v)
A = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
basis = calculate_interpolation_polynomials(p, A)
@test basis[1] == :(1.0 + -u + -v)
@test basis[2] == :(+u)
@test basis[3] == :(+v)

dbasis = calculate_interpolation_polynomial_derivatives(basis, 2)
@test isapprox(dbasis[1,1], -1.0)
@test isapprox(dbasis[2,1], -1.0)
@test isapprox(dbasis[1,2], 1.0)
@test isapprox(dbasis[2,2], 0.0)
@test isapprox(dbasis[1,3], 0.0)
@test isapprox(dbasis[2,3], 1.0)
