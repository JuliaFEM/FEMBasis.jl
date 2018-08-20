# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Test

@testset "FEMBasis.jl" begin
    @testset "test substitution" begin include("test_substitute.jl") end
    @testset "create Vandermonde matrix" begin include("test_vandermonde.jl") end
    @testset "create interpolation polynomials" begin include("test_generate_polynomials.jl") end
    @testset "test_create_basis" begin include("test_create_basis.jl") end
    @testset "test_lagrange_elements" begin include("test_lagrange_elements.jl") end
    @testset "test_nurbs_elements" begin include("test_nurbs_elements.jl") end
    @testset "test_math" begin include("test_math.jl") end
end
