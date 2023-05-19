#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

function rand_poly_dist(T, 
    vars::AbstractVector, 
    d::Integer,
    dist::String; 
    homogeneous::Bool = false
)
    """
    Create a random dense polynomial of degree `d` in the given variables `variables`.
    Each coefficient is sampled independently from a Normal(0, 1) or a Uniform(-0.5,0.5).
    """
    M = monomials(vars, d; affine = !homogeneous)
    if dist == "normal"
        sum(randn(T, length(M)) .* M)
    else dist == "uniform"
        sum((rand(T, length(M)).-0.5).*M)
    end
end

function randomsystem(vars, d, n, dist)
    """
    Build system of random polynomials
    """
    #initialize place holder for polynomial system
    equations = []
    #construct system
    for i in 1:n
        append!(equations, rand_poly_dist(Float64, vars, d, dist))
    end
    System(equations)
end

function matrix2poly(r, A...)
    """
    Map matrices to polynomials
    """
    #construct vector of variables
    n = length(r)
    @var x[1:n]
    #construct equation of species i
    for i in 1:n
        #get tensor of interactions for species i
        B_i = B[repeat([Colon()],d-1)...,i]
    end
end

function main(div_vec, hois_vec, n_sim, var, dist, stability=false)
    """
    Run bulk of simulations. Each simulation computes the number of real, and feasible
    equilibria, and their local stability, for a system of n species with interactions
    up to order d, when interaction coefficients are iid variables centered at 0, sampled
    from a certain distribution (so far, it can be gaussian or uniform).
    """
    n_hoi = length(hois_vec)
    n_div = length(div_vec)
    #set the number of columns depending on what is to be stored
    if stability n_col = 6 else n_col = 5 end
    #preallocate matrix to store number of zeros
    n_eq_mat = zeros(Float64, 1, n_col)
    #loop over higher order dimensions
    for d in hois_vec
        print("Dimension of HOIs: ")
        println(d)
        #loop over diversities
        for n in div_vec
            print("Diversity: ")
            println(n)
            #declare dynamic variables
            @var x[1:n]
            #solve n_sim realizations of communities with i spp
            for s in 1:n_sim
                system = randomsystem(x, d, n, dist)
                #show solving progress only if slow
                if n^d > 1024 show = true else show = false end
                #solve system and get real solutions
                real_sols = real_solutions(solve(system, show_progress = show))
                n_real = length(real_sols)
                #get positive solutions
                pos_sols = filter(s -> all(s .> 0), real_sols)
                n_pos = length(pos_sols)
                if stability && n_pos > 0 
                    #get largest eigenvalues of jacobian evaluated at all feasible equilibria
                    lambda_max_vec = zeros(Float64, n_pos)
                    for i in 1:n_pos
                        j_eval_i = jacobian(system, pos_sols[i])
                        lambda_max_vec[i] = maximum(real(eigvals(j_eval_i)))
                    end
                    #store results
                    add_rows = hcat(repeat([d n s n_real n_pos], n_pos), lambda_max_vec)
                    n_eq_mat = vcat(n_eq_mat, add_rows)
                elseif stability && n_pos == 0
                    #store a 1 flag for eigenvalue in the absence of feasible equilibria
                    add_rows = hcat([d n s n_real n_pos], 1)
                    n_eq_mat = vcat(n_eq_mat, add_rows)
                else 
                    #store results without stability measures
                    add_rows = [d n s n_real n_pos]
                    n_eq_mat = vcat(n_eq_mat, add_rows)
                end
            end
        end
    end
    return n_eq_mat
end

#set parameters
hoi_vec = [1 2 3 4 5] #polynomial degrees
div_vec = [2 3 4 5] #number of species
n_sim = 1000 #number of simulations
var = 1
dist = "normal"
stability = true
#run simulations
data = main(div_vec, hoi_vec, n_sim, var, dist, stability)
#save data
max_hoi = string(maximum(hoi_vec))
max_div = string(maximum(div_vec))
stab = string(stability)
output_name = "dim_"*max_hoi*"_div_"*max_hoi*"_s_"*string(n_sim)*"_"*dist*"_stab"*stab
writedlm("../data/expected_n_roots_"*output_name*".csv", data)