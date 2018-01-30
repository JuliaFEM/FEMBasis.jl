# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Base.Test
using FEMBasis
#using FEMBasis: calculate_basis_coefficients, calculate_interpolation_polynomials,
#                calculate_interpolation_polynomial_derivatives, create_lagrange_basis

@testset "FEMBasis.jl" begin

@testset "Calculate interpolation polynomial matrix" begin
    p = :(1 + u + v)
    X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    A = FEMBasis.calculate_basis_coefficients(p, X)
    A_expected = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
    @test isapprox(A, A_expected)
end

@testset "Calculate interpolation polynomials" begin
    p = :(1 + u + v)
    A = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
    equations = FEMBasis.calculate_interpolation_polynomials(p, A)
    @test equations[1] == :(1.0 + -u + -v)
    @test equations[2] == :(+u)
    @test equations[3] == :(+v)
end


@testset "Calculate derivatives of interpolation polynomials" begin
    basis = [:(1.0 - u - v), :(1.0u), :(1.0v)]
    dbasis = FEMBasis.calculate_interpolation_polynomial_derivatives(basis, 2)
    @test isapprox(dbasis[1][1], -1.0)
    @test isapprox(dbasis[1][2], -1.0)
    @test isapprox(dbasis[2][1], 1.0)
    @test isapprox(dbasis[2][2], 0.0)
    @test isapprox(dbasis[3][1], 0.0)
    @test isapprox(dbasis[3][2], 1.0)
end

@testset "test create basis, N and dN given" begin
    X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    basis = [:(1.0 - u - v), :(1.0u), :(1.0v)]
    dbasis = Vector[[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]
    code = FEMBasis.create_basis(:TestTriangle, "test triangle", X, basis, dbasis)
    eval(code)
    N = zeros(1,size(TestTriangle, 2))
    X = get_reference_element_coordinates(TestTriangle)
    for i=1:length(TestTriangle)
        eval_basis!(TestTriangle, N, X[i])
        N_expected = zeros(1, length(TestTriangle))
        N_expected[i] = 1.0
        @test isapprox(N, N_expected)
    end
    dN = zeros(size(TestTriangle)...)
    eval_dbasis!(TestTriangle, dN, X[1])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[2])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[3])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
end

@testset "test create basis, N given" begin
    X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    basis = ( :(1 - u - v), :(1.0u), :(1.0v) )
    code = FEMBasis.create_basis(:TestTriangle, "test triangle", X, basis)
    eval(code)
    N = zeros(1,size(TestTriangle, 2))
    X = get_reference_element_coordinates(TestTriangle)
    for i=1:length(TestTriangle)
        eval_basis!(TestTriangle, N, X[i])
        N_expected = zeros(1, length(TestTriangle))
        N_expected[i] = 1.0
        @test isapprox(N, N_expected)
    end
    dN = zeros(size(TestTriangle)...)
    eval_dbasis!(TestTriangle, dN, X[1])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[2])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[3])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
end

@testset "test create basis, ansatz polynomial given" begin
    X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    p1 = :(1 + u + v)
    p2 = "1 + u + v"
    code1 = FEMBasis.create_basis(:TestTriangle, "test triangle", X, p1)
    code2 = FEMBasis.create_basis(:TestTriangle, "test triangle", X, p2)
    @test code1 == code2
    eval(code1)
    N = zeros(1,size(TestTriangle, 2))
    X = get_reference_element_coordinates(TestTriangle)
    for i=1:length(TestTriangle)
        eval_basis!(TestTriangle, N, X[i])
        N_expected = zeros(1, length(TestTriangle))
        N_expected[i] = 1.0
        @test isapprox(N, N_expected)
    end
    dN = zeros(size(TestTriangle)...)
    eval_dbasis!(TestTriangle, dN, X[1])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[2])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
    eval_dbasis!(TestTriangle, dN, X[3])
    @test isapprox(dN, [-1.0 1.0 0.0; -1.0 0.0 1.0])
end

end
