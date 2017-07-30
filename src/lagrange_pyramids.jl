# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

#=
code = create_lagrange_basis(
    :Pyr5,
    "5 node linear pyramid element",
    "",
    (
     (-1.0, -1.0),
     ( 1.0, -1.0),
     ( 1.0,  1.0),
     (-1.0,  1.0)
    )
   )
eval(code)

"Pyr5" => "5 node linear pyramid element",


function get_reference_coordinates(::Type{Pyr5})
    Vector{Float64}[
        [-1.0,-1.0,-1.0], # N1
        [ 1.0,-1.0,-1.0], # N2
        [ 1.0, 1.0,-1.0], # N3
        [-1.0, 1.0,-1.0], # N4
        [ 0.0, 0.0, 1.0]] # N5
end

function get_interpolation_polynomial(::Type{Pyr5}, xi)
    [
        1.0/8.0*(1.0-1.0*xi[1])*(1.0-1.0*xi[2])*(1.0-1.0*xi[3])
        1.0/8.0*(1.0+1.0*xi[1])*(1.0-1.0*xi[2])*(1.0-1.0*xi[3])
        1.0/8.0*(1.0+1.0*xi[1])*(1.0+1.0*xi[2])*(1.0-1.0*xi[3])
        1.0/8.0*(1.0-1.0*xi[1])*(1.0+1.0*xi[2])*(1.0-1.0*xi[3])
        1.0/2.0*(1.0+xi[3])
    ]'
end

function get_interpolation_polynomial(::Type{Pyr5}, xi, ::Type{Val{:partial_derivatives}})
    [
        -0.125*(1.0-xi[2])*(1.0-xi[3])   0.125*(1.0-xi[2])*(1.0-xi[3])   0.125*(1.0+xi[2])*(1.0-xi[3])  -0.125*(1.0+xi[2])*(1.0-xi[3])  0.0
        -0.125*(1.0-xi[1])*(1.0-xi[3])  -0.125*(1.0+xi[1])*(1.0-xi[3])   0.125*(1.0+xi[1])*(1.0-xi[3])   0.125*(1.0-xi[1])*(1.0-xi[3])  0.0
        -0.125*(1.0-xi[1])*(1.0-xi[2])  -0.125*(1.0+xi[1])*(1.0-xi[2])  -0.125*(1.0+xi[1])*(1.0+xi[2])  -0.125*(1.0-xi[1])*(1.0+xi[2])  0.5
    ]
end
=#
