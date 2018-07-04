# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Calculus

function get_reference_element_coordinates end
function eval_basis! end
function eval_dbasis! end

function get_dim(p::Expr)
    pargs = p.args[2:end]
    :w in pargs && return 3
    :v in pargs && return 2
    :u in pargs && return 1
end

"""
    function calculate_basis_coefficients(polynomial::String, coordinates::Vararg{Tuple})

Calculate "interpolate coefficient matrix" for some polynomial p.

# Examples
That is, if we have polynomial p = 1 + u + v and coordinates (0,0), (1,0), (0,1),
we find A*p such that first row is the first coordinate, second row is second
coordinate and so on:
```julia
julia> p = "1 + u + w"
julia> X = ((0.0,0.0), (1.0,0.0), (0.0,1.0))
julia> calculate_basis_coefficient(p, X)
[1.0 0.0 0.0 # <-- p(0.0,0.0) = 1.0     = [1.0 0.0 0.0] * [1.0, u, v]
 1.0 1.0 0.0 # <-- p(1.0,0.0) = 1.0 + u = [1.0 1.0 0.0] * [1.0, u, v]
 1.0 0.0 1.0] # <- p(0.0,1.0) = 1.0 + v = [1.0 0.0 1.0] * [1.0, u, v]
```
"""
function calculate_basis_coefficients(p::Expr, X::Tuple)
    pargs = p.args[2:end]
    dim = get_dim(p)
    nbasis = length(X)
    sandbox = Module(:__SANDBOX__)
    vars = (:u, :v, :w)

    # calculate coefficient matrix A and inversion
    A = zeros(nbasis, nbasis)
    for i=1:nbasis
        for j=1:nbasis
            for (k,c) in enumerate(X[j])
                code = :( $(vars[k]) = $c )
                eval(sandbox, code)
            end
            A[j,i] = eval(sandbox, pargs[i])
        end
    end
    return A 
end

"""
"""
function calculate_interpolation_polynomials(p::Expr, A::Matrix)
    pargs = p.args[2:end]
    dim = get_dim(p)
    nbasis = size(A, 1)
    invA = inv(A)
    equations = Expr[]
    for j=1:nbasis
        equation = Expr(:call, :+)
        for i=1:nbasis
            isapprox(invA[i,j], 0.0) && continue
            if typeof(pargs[i]) <: Int && pargs[i] == 1
                push!(equation.args, invA[i,j])
            else
                if isapprox(invA[i,j], 1.0)
                    push!(equation.args, :($(pargs[i])))
                elseif isapprox(invA[i,j], -1.0)
                    push!(equation.args, :(-$(pargs[i])))
                else
                    push!(equation.args, :($(invA[i,j]) * $(pargs[i])))
                end
            end
        end
        push!(equations, equation)
    end
    return equations
end

"""
Calculate derivatives of basis functions with respect to parameters u, v, w.
"""
function calculate_interpolation_polynomial_derivatives(basis, dim)
    equations = Vector[]
    vars = [:u, :v, :w]
    nbasis = length(basis)
    for j=1:nbasis
        deq = []
        for k=1:dim
            #println("∂($(basis[j])) / ∂$(vars[k])")
            der = simplify(differentiate(basis[j], vars[k]))
            push!(deq, der)
        end
        push!(equations, deq)
    end
    return equations
end

function create_basis(name::Symbol, description::String, X::NTuple{nbasis, NTuple{dim, Float64}}, basis, dbasis) where {nbasis,dim}

    Q = Expr(:block)
    for i=1:nbasis
        push!(Q.args, :(N[$i] = $(basis[i])))
    end

    V = Expr(:block)
    for i=1:nbasis
        for k=1:dim
            push!(V.args, :(dN[$k,$i] = $(dbasis[i][k])))
        end
    end

    # FIXME: this can be done better somehow
    if dim == 1
        unpack = :((u,) = xi)
    elseif dim == 2
        unpack = :((u, v) = xi)
    else
        unpack = :((u, v, w) = xi)
    end

    code = quote
        mutable struct $name <: FEMBasis.AbstractBasis
        end

        function Base.size(::Type{$name})
            return ($dim, $nbasis)
        end

        function Base.size(::Type{$name}, j::Int)
            j == 1 && return $dim
            j == 2 && return $nbasis
        end

        function Base.length(::Type{$name})
            return $nbasis
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

function create_basis(name::Symbol, description::String, X::NTuple{nbasis, NTuple{dim, Float64}}, p::Expr) where {nbasis,dim}
    A = calculate_basis_coefficients(p, X)
    basis = calculate_interpolation_polynomials(p, A)
    dbasis = calculate_interpolation_polynomial_derivatives(basis, dim)
    return create_basis(name, description, X, basis, dbasis)
end

function create_basis(name::Symbol, description::String, X::NTuple{nbasis, NTuple{dim, Float64}}, p::String) where {nbasis,dim}
    create_basis(name, description, X, parse(p))
end

function create_basis(name::Symbol, description::String, X::NTuple{nbasis, NTuple{dim, Float64}}, basis_::NTuple{nbasis, String}) where {nbasis,dim}
    basis = [parse(b) for b in basis_]
    dbasis = calculate_interpolation_polynomial_derivatives(basis, dim)
    return create_basis(name, description, X, basis, dbasis)
end

function create_basis(name::Symbol, description::String, X::NTuple{nbasis, NTuple{dim, Float64}}, basis::NTuple{nbasis, Expr}) where {nbasis,dim}
    dbasis = calculate_interpolation_polynomial_derivatives(basis, dim)
    return create_basis(name, description, X, basis, dbasis)
end
