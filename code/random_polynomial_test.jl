#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles #to load and save files

function random_B(n)
    """
    Sample entries from a nx(n+1)x(n+1) tensor 
    from a gaussian distribution with mean 0 and sd 1
    """
    rng = MersenneTwister()
    randn(rng, Float64, (n+1,n+1, n))
end

function build_equation_i(x, B, i)
    """
    Build equation as in the paper
    """
    dot(x, B[:, :, i]*x)
end

function build_glv_i(x, r, A, B, i)
    """
    Build glv for species i
    """
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

function build_system(B, n)
    """
    Build system of polynomials
    """
    #declare dynamic variables
    @var x[1:n]
    #add a constant species
    #x = [x; 1]
    #initialize system of ODEs as empty list
    equations = []
    #construct set of dynamic equations
    for j in 1:n
        #eqn = build_equation_i(x, B, j)
        eqn = build_glv_i(x, r, A, B, i)
        append!(equations, eqn)
    end
    System(equations)
end

function count_feasible_roots(system)
    """
    Count number of positive solutions of polynomial system
    """
    #solve system
    res = solve(system, show_progress = false)
    #get only poisitive and real solutions
    valid_real_sols = filter(s -> all(s .> 0), real_solutions(res))
    #count them
    length(valid_real_sols)
end

function main(div_max, n_sim)
    """
    Get number of positive roots for n_sim simulations
    of n species, with n running from 3 to n_max
    """
    div_vec = 3:div_max
    n_div = length(div_vec)
    #preallocate storing matrix, initialize iterator
    n_eq_mat = zeros(n_div*n_sim, 3)
    it = 1
    #loop over diversities
    for i in 3:div_max
        #solve n_sim realizations of communities with i spp
        print("Diversity: ")
        println(i)
        for j in 1:n_sim
            B = random_B(i)
            system = build_system(B, i)
            n_zeros = count_feasible_roots(system)
            add_row = [i j n_zeros]
            n_eq_mat[it,:] = add_row
            it += 1
        end
    end
    n_eq_mat
end

div_max = 10
n_sim = 2000
data = main(div_max, n_sim)
writedlm("../data/expected_n_roots.csv", data)
