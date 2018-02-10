# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Base.Test
using FEMBasis

@testset "interpolate scalar field" begin
    # in unit square: T(X,t) = t*(1 + X[1] + 3*X[2] - 2*X[1]*X[2])
    B = Quad4()
    X = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    T = (1.0, 2.0, 3.0, 4.0)
    T_known(X) = 1 + X[1] + 3*X[2] - 2*X[1]*X[2]
    T_interpolated = interpolate(B, T, (0.0, 0.0))
    @test isapprox(T_interpolated, T_known((0.5, 0.5)))
end

@testset "calculate gradient dB/dX" begin
    B = Quad4()
    X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
    dBdX = grad(B, X, (0.0, 0.0))
    @test isapprox(dBdX, 1/2*[-1 1 1 -1; -1 -1 1 1])
end

@testset "interpolate gradient of scalar field" begin
    # in unit square: grad(T)(X) = [1-2X[2], 3-2*X[1]]
    B = Quad4()
    X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
    T = (1.0, 2.0, 3.0, 4.0)
    dTdX = grad(B, T, X, (0.0, 0.0))
    dTdX_expected(X) = [1-2*X[2] 3-2*X[1]]
    @test isapprox(dTdX, dTdX_expected([0.5, 0.5]))
end

@testset "interpolate vector field" begin
    # in unit square, u(X,t) = [1/4*t*X[1]*X[2], 0, 0]
    B = Quad4()
    u = Vector[[0.0, 0.0], [0.0, 0.0], [1/4, 0.0], [0.0, 0.0]]
    u_known(X) = [1/4*X[1]*X[2], 0]
    u_interpolated = interpolate(B, u, (0.0, 0.0))
    @test isapprox(u_interpolated, u_known((0.5, 0.5)))
end

@testset "interpolate gradient of vector field" begin
    # in unit square, u(X) = t*[X[1]*(X[2]+1), X[1]*(4*X[2]-1)]
    # => u_i,j = t*[X[2]+1 X[1]; 4*X[2]-1 4*X[1]]
    B = Quad4()
    X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])
    u = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])
    dudX = grad(B, u, X, (0.0, 0.0))
    dudX_expected(X) = [X[2]+1 X[1]; 4*X[2]-1 4*X[1]]
    @test isapprox(dudX, dudX_expected([0.5, 0.5]))
end

@testset "test BasisInfo" begin
    B = BasisInfo(Seg2)
    @test length(B) == 2
    @test size(B) == (1, 2)
    X = ((0.0,), (1.1))
    xi = (0.0,)
    eval_basis!(B, X, xi)
    @test isapprox(B.invJ, inv(B.J))
    @test isapprox(B.detJ, det(B.J))
    B = BasisInfo(Quad4)
    X = ((0.0,0.0), (1.1,0.0), (1.0,1.0), (0.0,1.0))
    xi = (0.0,0.0)
    eval_basis!(B, X, xi)
    @test isapprox(B.invJ, inv(B.J))
    @test isapprox(B.detJ, det(B.J))
    B = BasisInfo(Hex8)
    X = ((0.0,0.0,0.0), (1.1,0.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0),
         (0.0,0.0,1.0), (1.0,0.0,1.0), (1.0,1.0,1.0), (0.0,1.0,1.0))
    xi = (0.0,0.0,0.0)
    eval_basis!(B, X, xi)
    @test isapprox(B.invJ, inv(B.J))
    @test isapprox(B.detJ, det(B.J))
end

@testset "evaluate gradient using BasisInfo" begin
    B = BasisInfo(Quad4)
    X = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))
    xi = (0.0, 0.0)
    eval_basis!(B, X, xi)
    u = ((0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0))
    gradu = zeros(2, 2)
    grad!(B, gradu, u)
    gradu_expected = [1.5 0.5; 1.0  2.0]
    @test isapprox(gradu, gradu_expected)
end

@testset "test curves" begin
    bi = BasisInfo(Seg2)
    X1 = ((0.0,0.0,0.0), (1.0,0.0,0.0))
    X2 = ((0.0,0.0,0.0), (0.0,1.0,0.0))
    X3 = ((0.0,0.0,0.0), (0.0,0.0,1.0))
    xi = (0.0, 0.0)
    eval_basis!(bi, X1, xi)
    @test isapprox(bi.detJ, 0.5)
    eval_basis!(bi, X2, xi)
    @test isapprox(bi.detJ, 0.5)
    eval_basis!(bi, X3, xi)
    @test isapprox(bi.detJ, 0.5)
    @test isapprox(jacobian(Seg2(), X1, xi), [0.5 0.0 0.0])
    @test isapprox(jacobian(Seg2(), X2, xi), [0.0 0.5 0.0])
    @test isapprox(jacobian(Seg2(), X3, xi), [0.0 0.0 0.5])
end

@testset "test manifolds" begin
    bi = BasisInfo(Quad4)
    X1 = ((0.0,0.0,0.0), (1.0,0.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0))
    X2 = ((0.0,0.0,0.0), (0.0,1.0,0.0), (0.0,1.0,1.0), (0.0,0.0,1.0))
    X3 = ((0.0,0.0,0.0), (0.0,0.0,1.0), (1.0,0.0,1.0), (1.0,0.0,0.0))
    xi = (0.0, 0.0)
    eval_basis!(bi, X1, xi)
    @test isapprox(bi.detJ, 0.25)
    eval_basis!(bi, X2, xi)
    @test isapprox(bi.detJ, 0.25)
    eval_basis!(bi, X3, xi)
    @test isapprox(bi.detJ, 0.25)
    J1 = jacobian(Quad4(), X1, xi)
    J2 = jacobian(Quad4(), X2, xi)
    J3 = jacobian(Quad4(), X3, xi)
    @test isapprox(J1, [0.5 0.0 0.0; 0.0 0.5 0.0])
    @test isapprox(J2, [0.0 0.5 0.0; 0.0 0.0 0.5])
    @test isapprox(J3, [0.0 0.0 0.5; 0.5 0.0 0.0])
end
