# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using FEMBasis: subs
using Test

expression = :(1 + u + v + u*v^2)
data = (:u => 1.0, :v => 2.0)
result = subs(expression, data)
@debug "subs(expression, data) where" expression data
@test isapprox(result, 8.0)
