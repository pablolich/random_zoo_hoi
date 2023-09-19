using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix oducts
using DelimitedFiles #to load and save files

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
    variance(d::Int64, j_vec::Array)

calculate variance of coefficitent whose associated monomial has powers j_vec
"""
function variance(d::Int64, j_vec::Array)
    j_vec_expanded = [j_vec; (d - sum(j_vec))]
    multinomial(d, j_vec_expanded)
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
    buildstartsystem(n::Int64, d::Int64, vars::AbstractVector)

construct homotopy for all kss systems of n polynomials of degree d
"""
function buildstartsystem(D, vars::AbstractVector)
    System(vars .^ D .- 1, vars)
end

"""
    getstartsolutions(n::Int64, d::Int64)

construct all solutions for n equations of the form x^d-1 = 0
"""
function getstartsolutions(n::Int64, d::Int64)
    map(k -> [exp(im * 2 * pi * k[i] / d) for i in 1:n], Iterators.product(repeat([0:d-1],n)...))
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
    equations
end

"""
    buildall(n::Int64, d::Int64, nsim::Int64, vars::AbstractVector, rng::AbstractRNG)

construct nsim kss polynomial systems of n polynomials of degree d.
"""
function buildall(n::Int64, d::Int64, nsim::Int64, vars::AbstractVector, rng::AbstractRNG)
    allsystems = []
    #sample monomial structure
    M = monomials(vars, d)
    nmon = length(M)
    #loop over number of systems
    for i in 1:nsim
        system_i = buildsystem(M, nmon, vars, n, d, rng)
        if i == 1
            allsystems = system_i
        else
            allsystems = hcat(allsystems, system_i)
        end
    end
    return allsystems
end

"""
    is_feasible(r::PathResult)

return true if all components of r are positive
"""
function is_feasible(r::PathResult)
    return all(real(r.solution) .> 0) 
end

"""
    stopatfeasible(r::PathResult)

return true if r is feasible
"""
function stopatfeasible(r::PathResult)
    #check if solution  is real
    if is_real(r)
        #check if its feasible
        if is_feasible(r)
            return true
        else
            return false
        end
    else
        return false
    end
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

"""
    computefeasibility(systems, n::Int64, d::Int64, nsim::Int64)

determine the feasibility of all kss systems of a given n, d
"""
function computefeasibility(systems, n::Int64, d::Int64, nsim::Int64, vars::AbstractVector, save::Bool)
    results = []
    #loop over systems
    for i in 1:nsim
        #deal separately with the case of just one system
        if nsim == 1
            syst = systems
            powers, pars = 
            syst_fix = fix_parameters(syst)
        else
            syst = systems[:,i][1]
            syst_fix = fix_parameters(syst)
        end
        #solve system numerically
        println("solving system...(", i, ")")
        if parameter_homotopy
            #get one solution for the first iteration
            init_sol = []
            p1 = 
            if i == 1
                sols = solutions(solve(syst, 
                                       stop_early_cb = stopatfeasible,
                                       compile = false,
                                       start_system = :total_degree,
                                       threading = true,
                                       seed = UInt32(1),
                                       show_progress = false))
                init_sol = sols
            else

            end

        else
            sols = real_solutions(solve(syst, #track sols0 during the deformation of g to f
                                        stop_early_cb = stopatfeasible, #stop when a feasible solution is found
                                        compile = false, #not introduce compilation overhead
                                        start_system = :total_degree, #efficient way to start searching
                                        threading = true, #allow multithreading
                                        seed = UInt32(1), #seed for trackers
                                        show_progress = false))
        end
        #determine if there is a feasible solution
        pos_sols = filter(s -> all(s .> 0), sols)
        npos = length(pos_sols)
        #store
        append!(results, [n d npos])
    end
    results = reshape(results, 3, nsim)
    #save?
    if save
        open("../data/kss_simulations.csv", "a") do io
            writedlm(io, results', ' ')
        end
    end
end

"""
    getparameters(max_n::Int64, max_d::Int64, 
    comp_limit::Int64, specific_pairs::Bool)

construct tuple of parameter combinations for parameter sweep
"""
function getparameters(max_n, max_d, 
    comp_limit::Int64, specific_pairs::Bool)
    if specific_pairs
        n_d = hcat(max_n, max_d)
        nrow = size(n_d, 1)
        return [(n_d[i,1], n_d[i,2]) for i in 1:nrow]
    else
        return [(x, y) for x in 1:max_n, y in 1:max_d if y^x<comp_limit]
    end
end

"""
    parametersweep(nmax::Int64, dmax::Int64, nsim::Int64)

perform a parameter sweep with bounds given by nmax and dmax
"""
function parametersweep(nmax, dmax, nsim::Int64, seed::Int64)
    parameters = getparameters(nmax, dmax, 80000, true)
    n_pairs = length(parameters)
    #initialize random generator
    rng = MersenneTwister(seed)
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        println("System size ", n, " System degree ", d)
        #get variables
        @var x[1:n]
        #create all systems
        println("Building all systems...")
        systems = buildall(n, d, nsim, x, rng)
        computefeasibility(systems, n, d, nsim, x, true)
    end
end

seed = 1 

parametersweep([8, 7, 8],[6, 6, 5], 150, seed)