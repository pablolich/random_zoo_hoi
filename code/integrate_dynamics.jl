using SciMLBase
using LinearAlgebra
using HomotopyContinuation
using Random


function glv!(dy, y, s, W, t)
    Dy = Diagonal(y)
    dy =  Dy * (s - W*y)    
end

"""
    buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64, rng::AbstractRNG)

construct a kss system of n polynomials of degree d
"""
function buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64, rng::AbstractRNG)
    equations = []
    #loop over number of equations in each system
    for j in 1:n
        append!(equations, buildpoly(allmonomials, nmon, vars, d, rng))
    end
    System(equations)
end

"""
    variance(d::Int64, j_vec::Array)

calculate variance of coefficitent whose associated monomial has powers j_vec
"""
function variance(d::Int64, j_vec::Array)
    j_vec_expanded = [j_vec; (d - sum(j_vec))]
    multinomial(d, j_vec_expanded)
end

"""
    multinomialcoeff(n, kvec)

compute multinomial coefficient
"""
function multinomial(n::Int64, kvec::Array)
    num = factorial(n)
    den = prod([factorial(i) for i in kvec])
    num/den
end

"""
    buildpoly(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector)

given all the monomials, sample random coefficients for each and sum them to build a polynomial
"""
function buildpoly(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, d::Int64, rng::AbstractRNG)
    #initialize list of polynomial coefficients
    coefficient_list = zeros(Float16, nmon)
    for i in 1:nmon
        monomial_i = allmonomials[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, vars; expanded = true) #expanded reduces time
        #compute the variance of ith coefficient using mulitinomial coefficient
        var = variance(d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(var)*randn(rng, Float16)
    end
    sum(coefficient_list .* allmonomials)
end

"""
    storerow(i::Int64, storemat, tostore::Matrix{Int64})

storing results in a matrix
"""
function storerow(i::Int64, storemat, tostore::Matrix{Float64})
    if i == 1
        storemat = tostore
    else
        storemat = vcat(storemat, tostore)
    end
    return storemat
end

function getW(system::System, n::Int64, d::Int64, vars::AbstractVector)
    #initialize holders for coeffficients with first row
    O = []
    for i in 1:n
        #get equation i's exponents and coefficients
        exps, coeffs = exponents_coefficients(system.expressions[i], vars)
        #leave 0 exponents out
        coeffs0 = coeffs[1:(end-1)]
        #store
        O = storerow(i, O, reshape(coeffs0, (1, length(coeffs0))))
    end
    exp0 = exps[:, 1:(end-1)]
    println(exp0)
    W = transpose(exp0)*O
    return W
end

function main(vars::AbstractVector, n::Int64, d::Int64, rng::AbstractRNG)
    allmon = monomials(vars, d)
    #build system
    syst = buildsystem(allmon, length(allmon), vars, n, d, rng)
    #build matrix W
    W = getW(syst, n, d, vars)
end

seed = 1
rng = MersenneTwister(seed)
n = 2
d = 2
@var x[1:n]
#y0 = 