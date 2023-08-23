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
function buildpoly(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, d::Int64)
    #initialize list of polynomial coefficients
    coefficient_list = zeros(Float16, nmon)
    for i in 1:nmon
        monomial_i = allmonomials[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, vars; expanded = true) #expanded reduces time
        #compute the variance of ith coefficient using mulitinomial coefficient
        var = variance(d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(var)*randn(Float16)
    end
    sum(coefficient_list .* allmonomials)
end

"""
    buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64)

construct a kss system of n polynomials of degree d
"""
function buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64)
    equations = []
    #loop over number of equations in each system
    for j in 1:n
        append!(equations, buildpoly(allmonomials, nmon, vars, d))
    end
    System(equations)
end

"""
    buildall(n::Int64, d::Int64, nsim::Int64, vars::AbstractVector)

construct nsim kss polynomial systems of n polynomials of degree d.
"""
function buildall(n::Int64, d::Int64, nsim::Int64, vars::AbstractVector)
    allsystems = []
    #sample monomial structure
    M = monomials(vars, d)
    nmon = length(M)
    #loop over number of systems
    for i in 1:nsim
        system_i = buildsystem(M, nmon, vars, n, d)
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
    store(i::Int64, storemat, tostore::Matrix{Int64})

storing results in a matrix
"""
function storerow(i::Int64, storemat, tostore::Matrix{Int64})
    if i == 1
        storemat = tostore
    else
        storemat = vcat(storemat, tostore)
    end
    return storemat
end

"""
    computefeasibility(systems::Matrix{System}, n::Int64, d::Int64, nsim::Int64)

determine the feasibility of all kss systems of a given n, d
"""
function computefeasibility(systems::Matrix{System}, n::Int64, d::Int64, nsim::Int64, start::String)
    results = []
    #loop over systems
    for i in 1:nsim
        #solve system numerically
        syst = systems[:,i][1]
        if start == "total"
            sols = real_solutions(solve(syst, 
                                        stop_early_cb = stopatfeasible, #stop when a feasible solution is found
                                        compile = false, #not introduce compilation overhead
                                        start_system = :total_degree, #efficient way to start searching
                                        show_progress = false))
        else
            sols = real_solutions(solve(syst, 
                                        stop_early_cb = stopatfeasible, #stop when a feasible solution is found
                                        compile = false, #not introduce compilation overhead
                                        start_system = :polyhedral, #efficient way to start searching
                                        show_progress = false))
        end
        #determine if there is a feasible solution
        pos_sols = filter(s -> all(s .> 0), sols)
        npos = length(pos_sols)
        #store
        append!(results, [n d npos])
    end
    results = reshape(results, 3, nsim)
    #save
    open("../data/kss_simulations.csv", "a") do io
        writedlm(io, results', ' ')
    end
end

"""
    getparameters(max_n::Int64, max_d::Int64, 
    comp_limit::Int64, specific_pairs::Bool)

construct tuple of parameter combinations for parameter sweep
"""
function getparameters(max_n::Int64, max_d::Int64, 
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
function parametersweep(nmax::Int64, dmax::Int64, nsim::Int64, start::String)
    parameters = getparameters(nmax, dmax, 80000, true)
    n_pairs = length(parameters)
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        #get variables
        @var x[1:n]
        #create all systems
        systems = buildall(n, d, nsim, x)
        computefeasibility(systems, n, d, nsim, start)
    end 
end

@time parametersweep(6,6,5, "total")
@time parametersweep(6,6,5, "polyhedral")
#check if I can do the combinations [7, 6] and [8, 6]