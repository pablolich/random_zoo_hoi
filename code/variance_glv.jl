#build many systems of polynomials with n species and interaction degrees 1, 2, 3, 4

#GLV with three-way higher order interactions

using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using Statistics #to load and save files
using LinearAlgebra
using DelimitedFiles

function randomtensor(d, n, variance, dist)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    rng = MersenneTwister()
    if dist == "normal"
        return variance*randn(rng, Float64, repeat([n + 1], d)...)
    else dist == "uniform"
        return rand(rng, Float64, repeat([n + 1], d)...).-0.5
    end
end

function equation_i(T_i, vec)
    """
    Build GLV equation for species i with hois of dimension d  
    """
    #get dimension of T
    dim = ndims(T_i)
    n = size(T_i,1)
    #initialize total
    tot = 0
    if dim > 2
        #hois are more than three-way, build polynomial recursively
        for j in 1:n
            slices = repeat([Colon()], dim-1)
            tot += vec[j]*equation_i(T_i[slices...,j], vec)
        end
    elseif dim == 2
        #hois are three-way, compute quadratic form
        return dot(vec, T_i*vec)
    else
        #no hois, compute linear form
        return T_i'*vec
    end
    return tot
end

function getsystem(vars, B, d, n)
    """
    Build system of polynomials
    """
    #add a constant species
    vars = [vars; 1]
    #initialize system of ODEs as empty list
    equations = []
    #construct set of dynamic equations
    for i in 1:n
        #get tensor of interactions for species i
        B_i = B[repeat([Colon()],d-1)...,i]
        eqn = equation_i(B_i, vars)
        append!(equations, eqn)
    end
    System(equations)
end

function numbermonomials(d, n)
    """
    Number of monomials of degree d on the variables x1, ..., xn
    """
    return binomial(d+n-1, n-1)
end

function main(n, d, nsim, variance, dist)
    """
    Get variance of coefficients of polynomial form
    """
    @var x[1:n]
    println(x)
    #number of rows
    n_rows = sum([numbermonomials(i-1, n) for i in 1:d])
    #loop over simulations for each n
    coeffs_mat = zeros(Float64, n_rows, nsim)
    B = randomtensor(d, n, variance, dist)
    syst = getsystem(x, B, d, n)
    exp, coeffs = exponents_coefficients(syst.expressions[1], x)
    for i in 1:nsim
        B = randomtensor(d, n, variance, dist)
        syst = getsystem(x, B, d, n)
        #get coefficients
        coeffs_i = coefficients(syst.expressions[1], x)
        #store results
        coeffs_mat[:,i] = coeffs_i
        println(i)
    end
    mean_coeffs = mean(coeffs_mat, dims = 2)
    mean_square_coeffs = mean(coeffs_mat.^2, dims = 2)
    return(round.(mean_square_coeffs-mean_coeffs.^2), exp)
end

vars, exp = main(2, 5, 5000, 1, "normal")
writedlm("../data/variances.csv", hcat(vars, exp'))