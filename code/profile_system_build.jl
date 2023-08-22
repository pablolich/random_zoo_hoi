using Profile
using Random #to sample random tensors
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files
using HomotopyContinuation


function equation_i(T_i::Array, vec::Vector{Expression})
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

function randomtensor(d::Int64, n::Int64, rng::AbstractRNG)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    return randn(rng, Float64, repeat([n + 1], d)...)
end

function randomtensor(d::Int64, n::Int64, rng::AbstractRNG)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    return randn(rng, Float64, repeat([n + 1], d)...)
end

function multinomialcoeff(n, kvec)
    """
    Compute multinomial coefficient
    """
    num = factorial(n)
    den = prod([factorial(i) for i in kvec])
    return num/den
end

function variance_a(assumption::String, d::Integer, j_vec::Array)
    """
    Compute variance for specific assumption
    """
    if assumption=="kss"
        j_vec_expanded = [j_vec; (d - sum(j_vec))]
        variance = multinomialcoeff(d, j_vec_expanded)
    elseif assumption=="symmetric"
        variance = (multinomialcoeff(sum(j_vec), j_vec))^2
    else 
        print(assumption*" is not a valid assumption")
    end
    return variance
end

function getsystem(variables::AbstractVector{Variable}, 
    d::Int64, n::Int64, rng::AbstractRNG, fromtensor::Bool, 
    distribution::String)
    """
    Build system of polynomials
    """
    #initialize system for polynomial system
    equations = []
    println("Building system...")
    #build system from a random tensor
    if fromtensor
        #sample a random tensor
        B = randomtensor(d+1, n, rng)
        #add a constant species
        vars_homog = [variables; 1]
        #construct set of dynamic equations
        for i in 1:n
            #get tensor of interactions for species i
            B_i = B[repeat([Colon()],d)...,i]
            eqn = equation_i(B_i, vars_homog)
            append!(equations, eqn)
        end
    #build system from polynomial coefficients directly
    else
        for i in 1:n
            append!(equations, rand_poly_dist(Float16, variables, d, distribution, "kss"))
        end
    end
    System(equations)
end

function rand_poly_dist(T, 
    vars::AbstractVector, 
    d::Integer,
    distribution::String,
    assumption::String,
    homogeneous::Bool = false
)
    """
    Create a random dense polynomial of degree `d` in the given variables `variables`.
    Each coefficient is sampled independently from a Normal(0, var) or a Uniform(-0.5,0.5).
    """
    M = monomials(vars, d; affine = !homogeneous)
    print(M)
    n_terms = length(M)
    coefficient_list = zeros(Float16, n_terms)
    for i in 1:n_terms
        monomial_i = M[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, vars)
        if distribution == "gaussian"
            #compute the variance of ith coefficient using mulitinomial coefficient
            variance = variance_a(assumption, d, exponents)
            #sample ith coefficient from a gaussian with computed variance
            coefficient_list[i] = sqrt(variance)*randn(T)
        elseif distribution == "uniform"
            #compute limits of uniform ditribution
            b = uniformbounds(d, sum(exponents))
            #compute number of coefficients to be summed
            ncoeffs = multinomialcoeff(sum(exponents), exponents)
            #sample coefficients uniformly and sum them
            coefficient_list[i] = sum(rand(ncoeffs)*2*b.-b)
        else 
            print(distribution*" is not a valid distribution")
        end
    end
    sum(coefficient_list .* M)
end

function getparameters(max_n, max_d, 
    comp_limit::Int64, specific_pairs::Bool)
    """
    Creates parameter pairs while keeping complexity bounded
    """
    if specific_pairs
        n_d = hcat(max_n, max_d)
        nrow = size(n_d, 1)
        return [(n_d[i,1], n_d[i,2]) for i in 1:nrow]
    else
        return [(x, y) for x in 1:max_n, y in 1:max_d if y^x<comp_limit]
    end
end

function one_sweep()
    #set distribution from which to sample
    distribution = "gaussian"
    #get information about stability
    #form parameter pairs
    parameters = getparameters(5, 5, 80000, false)
    #set seed for parameter sweep
    rng = MersenneTwister(1)
    n_pairs = length(parameters)
    times = []
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        #declare dynamic variables
        @var x[1:n]
        #create and solve system
        t_recursive = @elapsed getsystem(x, d, n, rng, true, distribution) #fromtensor = true
        t_poly = @elapsed getsystem(x, d, n, rng, false, distribution) #fromtensor = false
        add_row = [n d t_recursive t_poly]
        if n_d == 1 
            #create object
            times = add_row
        else
            #append
            times = [times; add_row]
        end
    end
    return times
end

function many_sweeps(n_sweeps)
    many_times = []
    for i in 1:n_sweeps
        times = one_sweep()
        if i == 1 
            #create object
            many_times = times
        else
            #append
            many_times = [many_times; times]
        end
    end
    open("../data/times.csv", "w") do io
        writedlm(io, many_times, ' ')
    end
end

many_sweeps(2)