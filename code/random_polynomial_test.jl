#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles #to load and save files

function random_B(n_spp, dim, var)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #create tuple of tensor dimensions
    dims_t = tuple(repeat([n_spp], dim)...)
    rng = MersenneTwister()
    var*randn(rng, Float64, dims_t)
end

function build_equation_i(x, B, i)
    """
    Build equation as in the paper
    """
    dot(x, B[:, :, i]*x)
end

function build_eq_i(T_i, vec)
    #get dimension of T
    dim = ndims(T_i)
    n = size(T_i,1)
    #initialize total
    tot = 0
    if dim == 3
        return quadratic_form(T_i[:, :, i], vec)
    else 
        for j in 1:(n-1)
            tup_ind = tuple(repeat([Colon()], dim-1)...)
            tot += build_eq_i(T_i[tup_ind...,j], vec)
        end
    end
    return tot
end

function quadratic_form(A, x)
    dot(x, A*x)
end

function build_system(B, n)
    """
    Build system of polynomials
    """
    #declare dynamic variables
    @var x[1:n]
    #add a constant species
    x = [x; 1]
    #initialize system of ODEs as empty list
    equations = []
    #construct set of dynamic equations
    for j in 1:n
        eqn = build_equation_i(x, B, j)
        #eqn = build_glv_i(x, r, A, B, i)
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

function main(div_max, hois_d, n_sim, var)
    """
    Get number of positive roots for n_sim simulations
    of n species, with n running from 3 to n_max
    """
    div_vec = 3:div_max
    n_div = length(div_vec)
    #preallocate matrix to store number of zeros
    n_eq_mat = zeros(n_div*n_sim, 3)
    #initialize iterator
    it = 1
    #loop over diversities
    for i in 3:div_max
        #solve n_sim realizations of communities with i spp
        print("Diversity: ")
        println(i)
        for j in 1:n_sim
            B = random_B(i, var)
            system = build_system(B, i)
            n_zeros = count_feasible_roots(system)
            add_row = [i j n_zeros]
            n_eq_mat[it,:] = add_row
            it += 1
        end
    end
    n_eq_mat
end

div_max = 7
hoi_d = [2 3]
n_sim = 1000
var = 1
data = main(div_max, hoi_d, n_sim, var)
output_name = "dim_2"
writedlm("../data/expected_n_roots"*output_name*".csv", data)
