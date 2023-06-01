#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

function get_n_ds(max_n, max_d, comp_limit)
    return [(x, y) for x in 1:max_n, y in 1:max_d if y^x<5000]
end

function randomtensor(d, n, dist)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    rng = MersenneTwister()
    if dist == "normal"
        return randn(rng, Float64, repeat([n + 1], d)...)
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

function getsystem(vars, d, n, dist)
    """
    Build system of polynomials
    """
    #sample a random tensor
    B = randomtensor(d, n, dist)
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

###############################################################################

function multinomialcoeff(n, kvec)
    """
    Compute multinomial coefficient
    """
    num = factorial(n)
    den = prod([factorial(i) for i in kvec])
    return num/den
end

function variance_a(assumption, d, j_vec)
    """
    Compute variance for specific assumption
    """
    if assumption=="kss"
        j_vec_expanded = [j_vec; (d - sum(j_vec))]
        variance = multinomialcoeff(d, j_vec_expanded)
    elseif assumption=="symmetric"
        variance = (multinomialcoeff(sum(j_vec), j_vec))^2
    elseif assumption=="antisymmetric"
        n_perm = multinomialcoeff(sum(j_vec), j_vec)
        if iseven(n_perm)
            variance = 0
        else
            variance = 1
        end
    else 
        print(assumption*" is not a valid assumption")
    end
    return variance
end

function rand_poly_dist(T, 
    vars::AbstractVector, 
    d::Integer,
    assumption::String,
    homogeneous::Bool = false
)
    """
    Create a random dense polynomial of degree `d` in the given variables `variables`.
    Each coefficient is sampled independently from a Normal(0, var) or a Uniform(-0.5,0.5).
    """
    M = monomials(vars, d; affine = !homogeneous)
    n_terms = length(M)
    coefficient_list = zeros(Float64, n_terms)
    for i in 1:n_terms
        monomial_i = M[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, vars)
        #monomial_degree = degree(monomial_i)
        #compute the variance of ith coefficient using mulitinomial coefficient
        variance = variance_a(assumption, d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(variance)*randn(T)
    end
    #println(sum(coefficient_list .* M))
    sum(coefficient_list .* M)
end

function randomsystem(vars, d, n, assumption)
    """
    Build system of random polynomials
    """
    
    #initialize place holder for polynomial system
    equations = []
    #construct system
    for i in 1:n
        append!(equations, rand_poly_dist(Float64, vars, d, assumption))
    end
    System(equations)
end

function one_simulation(d, n, s, x, variance, dist, assumption)
    """
    Computes the number of real, and feasible equilibria, and (if asked) their local stability, 
    for a system of n species with interactions up to order d, 
    when interaction coefficients are iid variables centered at 0, sampled
    from a certain distribution (so far, it can be gaussian or uniform).
    """
    syst = randomsystem(x, d-1, n, assumption)
    #syst = getsystem(x, d, n, dist)
    #solve system and get real solutions
    real_sols = real_solutions(solve(syst, show_progress=false))
    nsol = length(real_sols)
    #get positive solutions
    pos_sols = filter(s -> all(s .> 0), real_sols)
    npos = length(pos_sols)
    if stability && npos > 0 
        #get largest eigenvalues of jacobian evaluated at all feasible equilibria
        lambda_max_vec = zeros(Float64, npos)
        for i in 1:npos
            j_eval_i = jacobian(syst, pos_sols[i])
            lambda_max_vec[i] = maximum(real(eigvals(j_eval_i)))
        end
        add_rows = hcat(repeat([d n s nsol npos], npos), lambda_max_vec)
    elseif stability && npos == 0
        #store a 1 flag for eigenvalue in the absence of feasible equilibria
        add_rows = hcat([n d nsol npos], 1)
    else 
        #store results without stability measures
        add_rows = [n d nsol npos]
    end
    return add_rows
end

function main(n_ds, n_sim, variance, dist, stability, assumption, merge, save_folder)
    """
    Run bulk of simulations
    """
    n_pairs = length(n_ds)
    #set the number of columns depending on what is to be stored
    if stability n_col = 5 else n_col = 4 end
    #loop over order of interactions and diversity pairs
    println("Running simulations under "*assumption*" case")
    for n_d in 1:n_pairs
        n = n_ds[n_d][1]
        d = n_ds[n_d][2]+1
        print("Dimension of HOIs: ")
        println(d)
        print("Diversity: ")
        println(n)
        #declare dynamic variables
        @var x[1:n]
        #preallocate matrix to store results
        n_eq_mat = zeros(Float64, n_sim, n_col)
        #loop over simulations for each n
        println("Simulation number:")
        for s in 1:n_sim
            #print progress
            if s==n_sim println(" ", s) else print(" ", s) end
            #store results
            add_rows = one_simulation(d, n, s, x, variance, dist, assumption)
            n_eq_mat[s,:] = add_rows
        end
        #after all simulations have ended, save a copy in the safe directory
        writedlm("../data/"*save_folder*"/n_"*string(n)*"_d_"*string(d)*".csv", n_eq_mat, ' ')
        #if merge is true, add a copy to the merging files situation
        if merge
            open("../data/"*save_folder*"/n_"*string(n)*"_d_"*string(d)*".csv", "a") do io
                writedlm(io, n_eq_mat, ' ')
            end
        end
    end
end

#set parameters
n_ds = get_n_ds(4, 6, 5000)
n_sim = 1000 #number of simulations
var = 1
dist = "normal"
stability = false
assumption = "symmetric"
save_folder = "kss_polynomials"
merge = false
#run simulations
@time main(n_ds, n_sim, var, dist, stability, assumption, merge, save_folder)
#save data
#max_deg = string(maximum(deg_vec))
#max_div = string(maximum(div_vec))
#stab = string(stability)
#output_name = "deg_"*max_deg*"_div_"*max_div*"_s_"*string(n_sim)*"_"*dist*"_stab"*stab
#writedlm("../data/expected_n_roots_"*output_name*".csv", data, delim="")