# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using FEMBasis
using FEMBasis: create_basis
using Test

X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
basis = [:(1.0 - u - v), :(1.0u), :(1.0v)]
dbasis = Vector[[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]

code = create_basis(:TestTriangle1, "test triangle 1", X, basis, dbasis)
eval(code)
N = zeros(1,size(TestTriangle1, 2))
X = get_reference_element_coordinates(TestTriangle1)
for i=1:length(TestTriangle1)
    eval_basis!(TestTriangle1, N, X[i])
    N_expected = zeros(1, length(TestTriangle1))
    N_expected[i] = 1.0
    @test isapprox(N, N_expected)
end
dN = zeros(size(TestTriangle1)...)
eval_dbasis!(TestTriangle1, dN, X[1])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle1, dN, X[2])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle1, dN, X[3])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])


code = create_basis(:TestTriangle2, "test triangle 2", X, basis)
eval(code)
N = zeros(1,size(TestTriangle2, 2))
X = get_reference_element_coordinates(TestTriangle2)
for i=1:length(TestTriangle2)
    eval_basis!(TestTriangle2, N, X[i])
    N_expected = zeros(1, length(TestTriangle2))
    N_expected[i] = 1.0
    @test isapprox(N, N_expected)
end
dN = zeros(size(TestTriangle2)...)
eval_dbasis!(TestTriangle2, dN, X[1])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle2, dN, X[2])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle2, dN, X[3])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])


p = :(1 + u + v)
code = create_basis(:TestTriangle3, "test triangle 3", X, p)
eval(code)
N = zeros(1, size(TestTriangle3, 2))
X = get_reference_element_coordinates(TestTriangle3)
for i=1:length(TestTriangle3)
    eval_basis!(TestTriangle3, N, X[i])
    N_expected = zeros(1, length(TestTriangle3))
    N_expected[i] = 1.0
    @test isapprox(N, N_expected)
end
dN = zeros(size(TestTriangle3)...)
eval_dbasis!(TestTriangle3, dN, X[1])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle3, dN, X[2])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
eval_dbasis!(TestTriangle3, dN, X[3])
@test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
