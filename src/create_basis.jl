# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Calculus

import Base.parse

function get_reference_element_coordinates end
function eval_basis! end
function eval_dbasis! end

function calculate_interpolation_polynomials(p, V)
    basis = []
    first(p.args) == :+ || error("Use only summation between terms of polynomial")
    args = p.args[2:end]
    N = size(V, 1)
    b = zeros(N)
    for i in 1:N
        fill!(b, 0.0)
        b[i] = 1.0
        solution = V \ b
        N = Expr(:call, :+)
        for (ai, bi) in zip(solution, args)
            isapprox(ai, 0.0) && continue
            push!(N.args, simplify( :( $ai * $bi ) ))
        end
        push!(basis, N)
    end
    return basis
end

function calculate_interpolation_polynomial_derivatives(basis, D)
    vars = [:u, :v, :w]
    dbasis = []
    for N in basis
        partial_derivatives = []
        for j in 1:D
            push!(partial_derivatives, simplify(differentiate(N, vars[j])))
        end
        push!(dbasis, partial_derivatives)
    end
    return dbasis
end

function create_basis(name, description, X::NTuple{N, NTuple{D, T}}, p::Expr) where {N, D, T}
    @debug "create basis given antsatz polynomial" name description X p
    V = vandermonde_matrix(p, X)
    basis = calculate_interpolation_polynomials(p, V)
    return create_basis(name, description, X, basis)
end

function create_basis(name, description, X::NTuple{N, NTuple{D, T}}, basis) where {N, D, T}
    @debug "create basis given basis functions" name description X basis
    dbasis = calculate_interpolation_polynomial_derivatives(basis, D)
    return create_basis(name, description, X, basis, dbasis)
end

function create_basis(name, description, X::NTuple{N, NTuple{D, T}}, basis, dbasis) where {N, D, T}

    @debug "create basis given basis functions and derivatives" name description X basis dbasis

    Q = Expr(:block)
    for i=1:N
        push!(Q.args, :(N[$i] = $(basis[i])))
    end

    V = Expr(:block)
    for i=1:N
        for k=1:D
            push!(V.args, :(dN[$k,$i] = $(dbasis[i][k])))
        end
    end

    if D == 1
        unpack = :((u,) = xi)
    elseif D == 2
        unpack = :((u, v) = xi)
    else
        unpack = :((u, v, w) = xi)
    end

    code = quote
        struct $name <: FEMBasis.AbstractBasis
        end

        function Base.size(::Type{$name})
            return ($D, $N)
        end

        function Base.size(::Type{$name}, j::Int)
            j == 1 && return $D
            j == 2 && return $N
        end

        function Base.length(::Type{$name})
            return $N
        end

        function FEMBasis.get_reference_element_coordinates(::Type{$name})
            return $X
        end

        function FEMBasis.eval_basis!(::Type{$name}, N::Matrix{T}, xi) where T
            $unpack
            $Q
            return N
        end

        function FEMBasis.eval_dbasis!(::Type{$name}, dN::Matrix{T}, xi) where T
            $unpack
            $V
            return dN
        end
    end

    return code
end

