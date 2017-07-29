# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Base.Test
using FEMBasis
using FEMBasis: calculate_basis_coefficients, calculate_interpolation_polynomials,
                calculate_interpolation_polynomial_derivatives, create_lagrange_basis

@testset "Calculate interpolation polynomial matrix" begin
    p = :(1 + u + v)
    X = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    A = calculate_basis_coefficients(p, X)
    A_expected = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
    @test isapprox(A, A_expected)
end

@testset "Calculate interpolation polynomials" begin
    p = :(1 + u + v)
    A = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
    equations = calculate_interpolation_polynomials(p, A)
    @test equations[1] == "N[1] = 1.0 - 1.0 * u - 1.0 * v"
end

@testset "Calculate derivatives of interpolation polynomials" begin
    p = :(1 + u + v)
    A = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0]
    equations = calculate_interpolation_polynomial_derivatives(p, A)
    @test equations[1] == "dN[1,1] = -1.0"
end

@testset "test create basis" begin
    p = "1 + u + v + u*v"
    X = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
    code = create_lagrange_basis(:MyQuad4, "4 node quadrangle element", p, X)
    println(code)
    #=
    @create_basis(
        :MyQuad4,
        "4 node quadrangle element",
        "1 + u + v + u*v",
        ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)))
    =#
    eval(code)

    N = zeros(1,size(MyQuad4, 2))
    X = get_reference_element_coordinates(MyQuad4)
    for i=1:length(MyQuad4)
        eval_basis!(MyQuad4, N, X[i])
        N_expected = zeros(1, length(MyQuad4))
        N_expected[i] = 1.0
        @test isapprox(N, N_expected)
    end

    dN = zeros(size(MyQuad4)...)
    eval_dbasis!(MyQuad4, dN, X[1])
    @test isapprox(dN, [-0.5 0.5 0.0 0.0; -0.5 0.0 0.0 0.5])
    eval_dbasis!(MyQuad4, dN, X[2])
    @test isapprox(dN, [-0.5 0.5 0.0 0.0; 0.0 -0.5 0.5 0.0])
    eval_dbasis!(MyQuad4, dN, X[3])
    @test isapprox(dN, [0.0 0.0 0.5 -0.5; 0.0 -0.5 0.5 0.0])
    eval_dbasis!(MyQuad4, dN, X[4])
    @test isapprox(dN, [0.0 0.0 0.5 -0.5; -0.5 0.0 0.0 0.5])
end
