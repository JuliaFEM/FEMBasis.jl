# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Base.Test

using FEMBasis

@testset "Quad4" begin
    N = zeros(1,size(Quad4, 2))
    X = get_reference_element_coordinates(Quad4)
    for i=1:length(Quad4)
        eval_basis!(Quad4, N, X[i])
        N_expected = zeros(1, length(Quad4))
        N_expected[i] = 1.0
        @test isapprox(N, N_expected)
    end
end
