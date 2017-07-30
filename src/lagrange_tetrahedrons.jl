# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

code = create_lagrange_basis(
    :Tet4,
    "4 node linear tetrahedral element",
    "1 + u + v + w",
    (
     (0.0, 0.0, 0.0), # N1
     (1.0, 0.0, 0.0), # N2
     (0.0, 1.0, 0.0), # N3
     (0.0, 0.0, 1.0), # N4
    )
   )
eval(code)

code = create_lagrange_basis(
    :Tet10,
    "10 node quadratic tetrahedral element",
    "1 + u + v + w + u*v + v*w + w*u + u^2 + v^2 + w^2",
    (
     (0.0, 0.0, 0.0), # N1
     (1.0, 0.0, 0.0), # N2
     (0.0, 1.0, 0.0), # N3
     (0.0, 0.0, 1.0), # N4
     (0.5, 0.0, 0.0), # N5
     (0.5, 0.5, 0.0), # N6
     (0.0, 0.5, 0.0), # N7
     (0.0, 0.0, 0.5), # N8
     (0.5, 0.0, 0.5), # N9
     (0.0, 0.5, 0.5), # N10
    )
   )
eval(code)
