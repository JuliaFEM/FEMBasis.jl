# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

code = create_lagrange_basis(
    :Seg2,
    "2 node linear segment/line element",
    "1 + u",
    ( (-1.0,), ( 1.0,) )
   )
#println(code)
eval(code)

code = create_lagrange_basis(
    :Seg3,
    "3 node quadratic segment/line element",
    "1 + u + u^2",
    ( (-1.0,), (1.0,), (0.0,))
   )
eval(code)
