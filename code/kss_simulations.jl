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
        vari = variance(d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(vari)*randn(rng, Float16)
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
    HomotopyContinuation.System(equations)
end

"""
    numberofcoefficients(d::Int64, n::Int64)

Compute number of monomials of a dense polynomial of degree d on n variables 
"""
function numberofcoefficients(d::Int64, n::Int64)
    ncoeffs = 0
    for i in 0:n
        ncoeffs += binomial(d+i-1, d-1)
    end
    return ncoeffs
end

"""
    samplecoefficients(poly::Expression, vars::AbstractVector, d::Int64, n::Int64)

return coefficients of one n-variate polynomial of degree d with the kss distribution
"""
function samplecoefficients(poly::Expression, vars::AbstractVector, d::Int64, n::Int64, times::Int64, rng::AbstractRNG)
    #calculate number of coefficients to sample
    ncoeffs = numberofcoefficients(d, n)
    coefficient_list = []
    #get matrix of exponents
    exponents, coeffs = exponents_coefficients(poly, vars)
    #loop over each monomial
    for i in 1:ncoeffs
        #calculate appropriate variance
        vari = variance(d, exponents[:, i])
        #sample coefficients
        append!(coefficient_list, sqrt(vari)*randn(rng, Float64, times))
    end
    return coefficient_list
end

"""
    buildparsystem(vars::AbstractVector, n::Int64, d::Int64)

Build a system of n dense multivariate polynomials of degree d on variables 
x1, ..., xn, where all the coefficients are parameters c_{ij}, where i in 1:n, and
j in 1:t, where t is calculated as numberofcoefficients(d, n)
"""
function buildparsystem(vars::AbstractVector, n::Int64, d::Int64)
    equations = []
    coefficients = []
    for i in 1:n
        polyi, coeffi = dense_poly(vars, d, coeff_name = Symbol(:c, "_", string(i)))
        append!(equations, polyi)
        append!(coefficients, coeffi)
    end
    return(System(equations, parameters = convert(Vector{Variable}, coefficients)))
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
    stopatreal(r::PathResult)

    return true if r is real
"""
function stopatreal(r::PathResult)
    #check if solution  is real
    if is_real(r)
        return true
    else
        return false
    end
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
            syst = System(systems)
        else
            syst = System(systems[:,i][1])
        end
        #solve system numerically
        println("solving system...(", i, ")")
        sols = real_solutions(solve(syst, #track sols0 during the deformation of g to f
                                    stop_early_cb = stopatfeasible, #stop when a feasible solution is found
                                    compile = false, #not introduce compilation overhead
                                    start_system = :total_degree, #efficient way to start searching
                                    threading = true, #allow multithreading
                                    seed = UInt32(1), #seed for trackers
                                    show_progress = false))
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
    getcoefficients(system, vars, d)

returns a vector with coefficients of all polynomials in the system
"""
function getcoefficients(system, vars, d)
    nexpr = length(system.expressions)
    coefficients = []
    for i in 1:nexpr
        append!(coefficients, coeffs_as_dense_poly(system.expressions[i], vars, d))
    end
    return coefficients
end

function solveandsave(seedx)
    rng = MersenneTwister(seedx)
    for n in 1:4
        for d in 1:4
            println("Diversity: ", n, " Interaction order: ", d)
            #set variables
            @var x[1:n]
            #sample monomial structure
            M = monomials(x, d)
            nmon = length(M)
            result = []
            #create a system
            syst = buildsystem(M, nmon, x, n, d, rng) #d-1 is the polynomial degree
            #solve system and get real solutions
            result = solve(syst,
                        compile = false, #not introduce compilation overhead
                        start_system = :total_degree, #efficient way to start searching
                        threading = true, #allow multithreading
                        seed = UInt32(1), #seed for trackers
                        show_progress = true)
            #save it
            solsfilename = "solutions_n_"*string(n)*"_d_"*string(d)
            write_solutions("../data/startsystems/"*solsfilename*".csv", solutions(result))
            #get coefficients of the system
            coeffs = getcoefficients(syst, x, d)
            coeffsfilename = "coeffs_n_"*string(n)*"_d_"*string(d)
            open("../data/startsystems/"*coeffsfilename*".csv", "a") do io
                writedlm(io, coeffs', ' ')
            end
            #save coefficients
        end
    end
end

function loadparameters(n, d)
    coeffsfilename = "coeffs_n_"*string(n)*"_d_"*string(d)
    return vec(readdlm("../data/startsystems/"*coeffsfilename*".csv"))
end

"""
    computefeasibility2(n::Int64, d::Int64, nsim::Int64)

determine the feasibility of all kss systems of a given n, d
"""
function computefeasibility2(n::Int64, d::Int64, nsim::Int64, vars::AbstractVector, rng::AbstractRNG, save::Bool)
    results = []
    #create generalized polynomial system
    syst = buildparsystem(vars, n, d)
    #load coefficients and precomputed solutions
    pinit = loadparameters(n, d)
    startsols = read_solutions("../data/startsystems/solutions_n_"*string(n)*"_d_"*string(d)*".csv")
    #generate all parameters at once
    all_pars = [samplecoefficients(syst.expressions[1], vars, d, n, n, rng) for _ in 1:nsim]
    #loop over systems
    for i in 1:nsim
        #   sample parameters of final system
        pfinal = samplecoefficients(syst.expressions[1], vars, d, n, n, rng)
        #solve system numerically
        println("solving system...(", i, ")")
        resul2 = solve(syst, startsols; #track startsols
                        start_parameters = pinit, 
                        target_parameters = pfinal, #perform parameter homotopy
                        #stop_early_cb = stopatfeasible, #stop when a feasible solution is found
                        compile = false, #not introduce compilation overhead
                        #start_system = :total_degree, #efficient way to start searching
                        threading = true, #allow multithreading
                        #seed = UInt32(1), #seed for trackers
                        show_progress = true)
        solutions(result2)
        #check if succeded
        foundroot = is_success(path_results(result)[1])
        if foundroot 
            sols = real_solutions(result)
        else
            #deal separately with the case of just one system
            if nsim == 1
                syst = System(systems)
            else
                syst = System(systems[:,i][1])
            end
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
    parameters = getparameters(nmax, dmax, 80000, false)
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
        computefeasibility(systems, n, d, nsim, x, false)
    end
end

"""
    parametersweep2(nmax::Int64, dmax::Int64, nsim::Int64)

perform a parameter sweep with bounds given by nmax and dmax
"""
function parametersweep2(nmax, dmax, nsim::Int64, seed::Int64)
    parameters = getparameters(nmax, dmax, 80000, false)
    n_pairs = length(parameters)
    #initialize random generator
    rng = MersenneTwister(seed)
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        println("System size ", n, " System degree ", d)
        #get variables
        @var x[1:n]
        computefeasibility2(n, d, nsim, x, rng, false)
    end
end


seedx = 1
solveandsave(seedx)
seedx = 2
@time parametersweep(4, 4, 10, seedx)
@time parametersweep2(4, 4, 10, seedx)