# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Base.Test

using FEMBasis

function test_elements(elements::Symbol...)
    for i in 1:length(elements)
        element = getfield(FEMBasis, elements[i])
        N = zeros(1,size(element, 2))
        X = get_reference_element_coordinates(element)
        for i=1:length(element)
            eval_basis!(element, N, X[i])
            N_expected = zeros(1, length(element))
            N_expected[i] = 1.0
            @test isapprox(N, N_expected)
        end
    end
end

@testset "Continuous Lagrange elements" begin
    @testset "Segments" begin test_elements(:Seg2, :Seg3) end
    @testset "Quadrangles" begin test_elements(:Quad4, :Quad8, :Quad9) end
    @testset "Triangles" begin test_elements(:Tri3, :Tri6, :Tri7) end
    @testset "Tetrahedrons" begin test_elements(:Tet4, :Tet10) end
    @testset "Hexahedrons" begin test_elements(:Hex8, :Hex20, :Hex27) end
    @testset "Wedges" begin test_elements(:Wedge6) end
    # fixme: partition of unity is not satisfied
    #@testset "Pyramids" begin test_elements(:Pyr5, ) end
    @testset "Pyramids" begin test_elements(:Pyr5CA, ) end
end
