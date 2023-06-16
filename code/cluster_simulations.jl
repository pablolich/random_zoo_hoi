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

function getsystem(variables::AbstractVector{Variable}, 
    d::Int64, n::Int64, rng::AbstractRNG)
    """
    Build system of polynomials
    """
    #sample a random tensor
    B = randomtensor(d+1, n, rng)
    #add a constant species
    vars_homog = [variables; 1]
    #initialize system of ODEs as empty list
    equations = []
    #construct set of dynamic equations
    for i in 1:n
        #get tensor of interactions for species i
        B_i = B[repeat([Colon()],d)...,i]
        eqn = equation_i(B_i, vars_homog)
        append!(equations, eqn)
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

function parameter_sweep(parameters, rng::AbstractRNG, save_rows)
    """
    Perform one sweep over all parameter
    """
    n_pairs = length(parameters)
    #initialize matrix of results
    sweepresult=[]
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        #declare dynamic variables
        @var x[1:n]
        #create and solve system
        syst = getsystem(x, d, n, rng)
        #solve system and get real solutions
        real_sols = real_solutions(solve(syst, show_progress=true))
        #get number of real and positive solutions
        nsol = length(real_sols)
        pos_sols = filter(s -> all(s .> 0), real_sols)
        npos = length(pos_sols)
        add_rows = [d n nsol npos]
        #save output to local file
        if save_rows
            open("../data/parameter_sweeps.csv", "a") do io
                writedlm(io, add_rows, ' ')
            end
        end
        if n_d == 1 
            #create object
            sweepresult = add_rows
        else
            #append
            sweepresult = [sweepresult; add_rows]
        end
    end
    return sweepresult
end

function manysweeps(n_sweeps::Int64, seed::Int64)
    """
    Perform multiple parameter sweeps
    """
    parameters = getparameters(8, 6, 20000, false)
    #set seed for parameter sweep
    rng = MersenneTwister(seed)
    #initialize matrix for storing results
    manysweeps_result = []
    println("Simulation number:")  
    for sweep in 1:n_sweeps
        if sweep==n_sweeps println(" ", sweep) elseif rem(sweep, 1)==0 print(" ", sweep)  else end
        sweepresult = parameter_sweep(parameters, rng, false)
        if sweep == 1 
            manysweeps_result = sweepresult 
        else
            manysweeps_result = [manysweeps_result; sweepresult] 
        end
    end
    return manysweeps_result
end

#manysweeps(1, 1) #n_sweeps, seed  
@time manysweeps(parse(Int, ARGS[1]), parse(Int, ARGS[2]))