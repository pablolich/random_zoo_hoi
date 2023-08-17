using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

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

function uniformbounds(d, monomial_degree)
    return sqrt(3/factorial(d-monomial_degree))
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

function getsystem(variables::AbstractVector{Variable}, 
    d::Int64, n::Int64, rng::AbstractRNG, fromtensor::Bool, 
    distribution::String)
    """
    Build system of polynomials
    """
    #initialize system for polynomial system
    equations = []
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
            append!(equations, rand_poly_dist(Float16, vars, d, distribution, assumption))
        end
    end
    System(equations)
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

function is_feasible(r::PathResult)
    return all(real(r.solution) .> 0) 
end

function is_stable(r::PathResult, jac_mat::Matrix{Expression})
    #define variables here
    @var x[1:size(jac_mat,1)]
    #evaluate jacobian at equilibirum
    jac_eval = evaluate(jac_mat, x => r.solution)
    #get eigenvalues
    lambda_max = maximum(real(eigvals(jac_eval)))
    if lambda_max < 0
        return true
    else
        return false
    end
end

function findstable(r::PathResult)
    #check if solution  is real
    if is_real(r)
        #check if its feasible
        if is_feasible(r)
            #check if its stable
            if is_stable(r, jac_mat)
                return true
            else
                return false
            end
        else
            return false
        end
    else
        return false
    end
end

"""
    getallmodifications(n, k)

generate all posible vectors of n choose k elements
"""
function getallmodifications(n, k)
    return collect(combinations(1:n, k))
end

"""
    modifysystem(system, modification, variables)

given a system, a modified system is returned with the power of modification equations, increased by 1, 
by multiplying by the species in the vector modification
"""
function modifysystem(system, modification, variables)
    modsyst = copy(system)
    for i in modification
        modsyst.expressions[i] = variables[i]*modsyst.expressions[i]
    end    
    return modsyst
end

"""
    is_any_stable(result::Result, jacobian::Matrix{Expression})

given the solutions of a system, determine if any of them is stable
"""
function is_any_stable(result::Result, jacobian::Matrix{Expression})
    #get number of solutions
    nres = nresults(result)
    for i in 1:nres
        sol_stab = is_stable(result.path_results[i], jacobian)
        if sol_stab
            return true
        end
    end
    return false
end

"""
    findlargeststable(syst)

get size of largest stable sub-equilibira
"""
function findlargeststable(syst, jacobian)
    nspp = length(syst)
    x = syst.variables
    #loop over how many species we are allowed zero out
    for i in 1:nspp
        #find all the combinations nspp choose i
        vec_mods = getallmodifications(nspp, i)
        #sequentially modify systems in rows of vec_mod
        for modification in vec_mods
            #modify system 
            mod_sys = modifysystem(syst, modification, x)
            #solve the system
            res = solve(mod_sys)
            #check if any of solutions are stable
            if is_any_stable(res, jacobian)
                #return size of the equilibrium that is stable
                return  nspp - length(modification)
            end
        end
    end
    #if none found, return 0
    return 0
end

function parameter_sweep(parameters, rng::AbstractRNG, save_rows::Bool, 
    distribution::String, stability::Bool) #feed directly a distribution?
    """
    Perform one sweep over all parameter
    """
    n_pairs = length(parameters)
    n_max = maximum([parameters[i][1] for i in 1:n_pairs])
    #initialize matrices of results
    sweep_sum_stat = []
    sweep_largest_stab = []
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        #declare dynamic variables
        @var x[1:n]
        #create and solve system
        syst = getsystem(x, d, n, rng, true, distribution) #fromtensor = true
        #construct global variable, jacobian
        global jac_mat = jacobian(syst)
        result = solve(syst, stop_early_cb = findstable)
        #get real solutions
        real_sols = real_solutions(result)
        #get number of real, positive, and stable solutions
        nsol = length(real_sols)
        pos_sols = filter(s -> all(s .> 0), real_sols)
        npos = length(pos_sols)
        nstab = 0
        #there are no feasible equilibria
        if stability && npos == 0
            #no feasible equilibria implies 0 stable equilibria
            sum_stat = [n d nsol npos 0]
            #get largest feasible stable
            nstab = findlargeststable(syst, jac_mat)
            largest_stab = [n d nstab]
        #there are some fasible equilibria
        else stability && npos > 0 
            #check if any of the feasible are stable
            if is_any_stable(result, jac_mat)
                sum_stat = [n d nsol npos 1]
                largest_stab = [n d n]
            #check for the largest feasible stable sub-equilibria
            else
                sum_stat = [n d nsol npos 0]
                nstab = findlargeststable(syst, jac_mat)
                largest_stab = [n d nstab]
            end
        end
        #save output to local file
        if save_rows
            open("../data/parameter_sweeps_summary.csv", "a") do io
                writedlm(io, sum_stat, ' ')
            end
            open("../data/parameter_sweeps_largest_stab.csv", "a") do io
                writedlm(io, largest_stab, ' ')
            end
        end
        if n_d == 1 
            #create object
            sweep_sum_stat = sum_stat
            sweep_largest_stab = largest_stab
        else
            #append
            sweep_sum_stat = [sweep_sum_stat; sum_stat]
            sweep_largest_stab = [sweep_largest_stab; largest_stab]
        end
    end
    return sweep_sum_stat, sweep_largest_stab
end

function manysweeps(n_sweeps::Int64, seed::Int64)
    """
    Perform multiple parameter sweeps
    """
    #set distribution from which to sample
    distribution = "normal"
    #get information about stability
    #form parameter pairs
    parameters = getparameters(8, 6, 20000, false)
    #set seed for parameter sweep
    rng = MersenneTwister(seed)
    #initialize matrix for storing results
    manysweeps_summary = []
    manysweeps_largest_stab = []
    println("Simulation number:")  
    for sweep in 1:n_sweeps
        if sweep==n_sweeps println(" ", sweep) elseif rem(sweep, 1)==0 print(" ", sweep)  else end
        sweep_summary, sweep_largest_stab = parameter_sweep(parameters, rng, true, distribution, true)
        if sweep == 1 
            manysweeps_summary = sweep_summary 
            manysweeps_largest_stab = sweep_largest_stab
        else
            manysweeps_summary = [manysweeps_summary; sweep_summary] 
            manysweeps_largest_stab = [manysweeps_largest_stab; sweep_largest_stab]
        end
    end
    open("../data/sweeps/n_sweeps_"*string(n_sweeps)*"_seed_"*string(seed)*"_summary.csv", "w") do io
        writedlm(io, manysweeps_summary, ' ')
    end
    open("../data/sweeps/n_sweeps_"*string(n_sweeps)*"_seed_"*string(seed)*"_largest_stab.csv", "w") do io
        writedlm(io, manysweeps_largest_stab, ' ')
    end
    return manysweeps_result
end

manysweeps(2000, 1) #n_sweeps, seed  
#@time manysweeps(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
    