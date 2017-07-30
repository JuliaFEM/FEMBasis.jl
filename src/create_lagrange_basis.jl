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
    debug("Polynomial: $pargs")
    dim = get_dim(p)
    nbasis = length(X)
    debug("dim = $dim, nbasis = $nbasis")
    sandbox = Module(:__SANDBOX__)
    vars = (:u, :v, :w)

    # calculate coefficient matrix A and inversion
    A = zeros(nbasis, nbasis)
    for i=1:nbasis
        for j=1:nbasis
            for (k,c) in enumerate(X[j])
                code = :( $(vars[k]) = $c )
                debug("evaluate code $code")
                eval(sandbox, code)
            end
            #eval(sandbox, parse("u, v, w = (0.0, 0.0, 0.0)"))
            # u, v = coordinates[j]
            term = pargs[i]
            debug("evaluate term $term in sandbox")
            A[j,i] = eval(sandbox, term)
        end
    end
    debug("A ready: ")
    debug(A)
    return A 
end

"""
"""
function calculate_interpolation_polynomials(p::Expr, A::Matrix)
    pargs = p.args[2:end]
    dim = get_dim(p)
    nbasis = size(A, 1)
    invA = inv(A)
    # evaluate polynomials
    equations = String[]
    for j=1:nbasis
        terms = String[]
        for i=1:nbasis
            isapprox(invA[i,j], 0.0) && continue
            if typeof(pargs[i]) <: Int && pargs[i] == 1
                push!(terms, "$(invA[i,j])")
            else
                push!(terms, "$(invA[i,j]) * $(pargs[i])")
            end
        end
        eq = "N[$j] = " * join(terms, " + ")
        eq = replace(eq, " + -", " - ")
        push!(equations, eq)
    end
    debug(join(equations, "\n"))
    return equations
end

"""
Calculate derivatives of basis functions with respect to parameters u, v, w.
"""
function calculate_interpolation_polynomial_derivatives(p::Expr, A::Matrix)
    pargs = p.args[2:end]
    dim = get_dim(p)
    nbasis = size(A, 1)
    invA = inv(A)
    # evaluate polynomials
    equations = String[]
    vars = [:u, :v, :w]
    for j=1:nbasis
        for k=1:dim
            terms = String[]
            for i=1:nbasis
                isapprox(invA[i,j], 0.0) && continue
                dterm = differentiate(pargs[i], vars[k])
                if typeof(dterm) <: Int
                    if dterm == 0
                        continue
                    end
                    if dterm == 1
                        push!(terms, "$(invA[i,j])")
                    end
                else
                    push!(terms, "$(invA[i,j]) * $dterm")
                end
            end
            if length(terms) == 0
                push!(terms, "0.0")
            end
            eq = "dN[$k,$j] = " * join(terms, " + ")
            #eq = replace(eq, "1.0 * 1 ", "1.0 ")
            #eq = replace(eq, " + -", " - ")
            push!(equations, eq)
        end
    end
    debug(join(equations, "\n"))
    return equations
end

function create_lagrange_basis(name::Symbol, description::String, ps::String, X::Tuple)

    p = parse(ps)
    dim = get_dim(p)
    nbasis = length(X)
    A = calculate_basis_coefficients(p, X)
    N = calculate_interpolation_polynomials(p, A)
    dN = calculate_interpolation_polynomial_derivatives(p, A)

    Q = Expr(:block)
    for i=1:nbasis
        push!(Q.args, parse(N[i]))
    end

    dN = reshape(dN, dim, nbasis)
    V = Expr(:block)
    for i=1:nbasis
        for k=1:dim
            push!(V.args, parse(dN[k, i]))
        end
    end

    # FIXME: this can be done better somehow
    if dim == 1
        xitype = :(Tuple{T})
        unpack = :((u,) = xi)
    elseif dim == 2
        xitype = :(Tuple{T,T})
        unpack = :((u, v) = xi)
    else
        xitype = :(Tuple{T,T,T})
        unpack = :((u, v, w) = xi)
    end

    code = quote
        type $name
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

        function FEMBasis.eval_basis!{T}(::Type{$name}, N::Matrix{T}, xi::$xitype)
            $unpack
            $Q
            return N
        end

        function FEMBasis.eval_dbasis!{T}(::Type{$name}, dN::Matrix{T}, xi::$xitype)
            $unpack
            $V
            return dN
        end
    end

    return code
end
