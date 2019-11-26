# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using FEMBasis
using Test

basis = NSeg()
@test size(basis) == (1, 2)
@test length(basis) == 2
N = zeros(2)
eval_basis!(basis, N, Vec(0.0, ))
@test isapprox(N, [0.5, 0.5])

basis = NSurf()
@test size(basis) == (2, 4)
@test length(basis) == 4
N = zeros(4)
eval_basis!(basis, N, Vec(0.0, 0.0))
@test isapprox(N, [0.25, 0.25, 0.25, 0.25])

basis = NSolid()
@test size(basis) == (3, 8)
@test length(basis) == 8
N = zeros(8)
eval_basis!(basis, N, Vec(0.0, 0.0, 0.0))
@test isapprox(N, [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125])
