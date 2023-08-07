using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files
using Symbolics #perform analytical calculations
using Roots

function evaluateglv(x::Array, A, r::Array)
    """
    function to evaluate original model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        response = []
        for j in 1:nspp
            if j == i
                append!(response, 0)
            else 
                append!(response, A[i,j]*x[j]/(1+x[j]))
            end
        end
        append!(res, r[i] - A[i,i]*x[i] - sum(response))
        println(res)
    end
    return res
end

function buildglvpoly(x::Array, A::Array, r::Array)
    """
    function to build polynomial form of the model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        terms = []
        factors = 1 + x
        factors[i] = 1
        for j in 1:nspp
            if j != i
                append!(terms, A[i,j]*x[j]*prod(factors)/(1+x[j]))
            end
            #try/catch to handle the particularly annoying case of one species
            try
                empty = !any(terms)
                if empty
                    append!(terms, 0)
                end
            catch
            end
        end
        append!(res, prod(factors)*(r[i] - A[i,i]*x[i]) - sum(terms))
    end
    return res
end

function certifysolutions(model, solutions, A, r)
    """
    Check if solutions are solutions to model or not
    """
    nsol = size(solutions, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(solutions[i], A, r)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        end
    end
    return cert_vec
end

"""
    interactionstructure(n, n_neighbor)

which species does species one interact with
"""
function interactionstructure(n, n_neighbor)
    first_row = []
    if n_neighbor == 0
        return first_row
    else
        #create vector with indices
        indices = collect(2:n)
        #reverse it
        indices_rev = reverse(indices)
        i = 1
        while i <= n_neighbor
            push!(first_row, indices[i])
            if length(first_row) == n_neighbor
                return first_row
            elseif length(first_row) > n_neighbor
                #how much longer is it
                longer = length(first_row) - n_neighbor
                for j in 0:(longer-1)
                    popat!(first_row, length(first_row)-j) 
                end
                return first_row
            else
                push!(first_row, indices_rev[i])
                i+=1
            end
        end
    end
    return first_row
end

"""
    getinteractionmatrix(interaction_strength, n, n_neighbor)

Build interaction matrix given the interaction structure
"""
function getinteractionmatrix(n, n_neighbor, interaction_strength) 
    """
    Build matrix of interactions
    """
    A = []
    #build first row of matrix
    ind_row = interactionstructure(n, n_neighbor)
    first_row = zeros(Int64, n)
    for i in ind_row
        first_row[i] = interaction_strength
    end
    #build circulant matrix
    for j in 1:n
        if j == 1
            A = first_row
        else
            A = hcat(A, circshift(first_row, j-1))
        end
    end
    return A' #transpose it
end

function buildfirstrow(n, n_neighbor, interacton_strength)
    #get non-zero indices
    ind_row = interactionstructure(n, n_neighbor)
    first_row = zeros(Int64, n)
    #set those positions to desired strength
    for i in ind_row
        first_row[i] = interacton_strength
    end
    return first_row
end

function buildcirculant(first_row)
    C = []
    n = length(first_row)
    for i in 1:n
        if j == 1
            C = first_row
        else
            C = hcat(A, circshift(first_row, j-1))
        end
    end
    return C' #transpose it
end

function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end

function getrepeatedeigs(eigenvalues)
    #evaluate eigenvalues numerically
    numeric_eigenvals = Vector{Float64}([])
    n_eigenvals = length(eigenvalues)
    for i in 1:n_eigenvals
        push!(numeric_eigenvals, Float64(real(substitute(eigs[i], Dict([a=>pi]))).val))
    end
    rounded_eigs = round.(numeric_eigenvals, digits = 3)
    rep_eigenvalues = nonunique(rounded_eigs)
    indices = findall(x->in(x, rep_eigenvalues), unique(rounded_eigs))
    return eigenvalues[indices]
end

function symbtofloat(symbolicexpression, symbolicvariable, numericvalue)
    return substitute(symbolicexpression, Dict([symbolicvariable=>numericvalue])).val
end
function findcritical(lambda, variable)
    findzero()
    symbtofloat(lambda, variable, )
end

#functions dealing with analytical calculations

function symmetricequilibrium(n_neighbor)
    @variables a
    return 1/2*(-n_neighbor*a + sqrt((n_neighbor*a)^2 + 4))
end

function jacobianelement(n_neighbor, bifurcation_par)
    return -4*bifurcation_par/(2-bifurcation_par*n_neighbor + sqrt(4+(bifurcation_par*n_neighbor)^2))^2
end

function rootsofunity(n)
    return exp(2*pi*im/n)
end

function jaceigenvals(c, n, n_neighbor)
    eigs = []
    sumind = interactionstructure(n, n_neighbor).-1
    for i in 1:n
        roots = 0
        for j in sumind
            roots += (rootsofunity(n))^((i-1)*j)
        end
        push!(eigs, -1 + c*roots)
    end
    return eigs
end