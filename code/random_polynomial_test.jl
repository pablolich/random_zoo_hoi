#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles #to load and save files

function random_B(dim, n_spp, var)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #create tuple of tensor dimensions
    dims_t = tuple(repeat([n_spp+1], dim)...)
    rng = MersenneTwister()
    var*randn(rng, Float64, dims_t)
end

function build_equation_i(x, B, i)
    """
    Build equation as in the paper
    """
    dot(x, B[:, :, i]*x)
end

function build_equation_i(T_i, i, vec)
    #get dimension of T
    dim = ndims(T_i)
    n = size(T_i,1)
    #initialize total
    tot = 0
    if dim == 3
        return quadratic_form(T_i[:,:,i], vec[1:n])
    else 
        for j in 1:(n-1)
            tup_ind = tuple(repeat([Colon()], dim-1)...)
            tot += build_equation_i(T_i[tup_ind...,j], i, vec)
        end
    end
    return tot
end

function quadratic_form(A, x)
    dot(x, A*x)
end

function build_system(B, d, n)
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
    for i in 1:n
        B_i = B[repeat([Colon()],d-1)...,i]
        eqn = build_equation_i(B_i, i, x)
        append!(equations, eqn)
    end
    System(equations)
end

function count_feasible_roots(system)
    """
    Count number of positive solutions of polynomial system
    """
    #solve system
    res = solve(system, show_progress = true)
    #get only poisitive and real solutions
    valid_real_sols = filter(s -> all(s .> 0), real_solutions(res))
    #count them
    length(valid_real_sols)
end

function main(div_vec, hois_vec, n_sim, var)
    """
    Get number of positive roots for n_sim simulations
    of n species, with n running from 3 to n_max
    """
    n_hoi = length(hois_vec)
    n_div = length(div_vec)
    #preallocate matrix to store number of zeros
    n_eq_mat = zeros(n_hoi*n_div*n_sim, 4)
    #initialize iterator
    it = 1
    #loop over higher order dimensions
    for d in hois_vec
        print("Dimension of HOIs: ")
        println(d)
        #loop over diversities
        for n in div_vec
            #solve n_sim realizations of communities with i spp
            print("Diversity: ")
            println(n)
            for s in 1:n_sim
                B = random_B(d, n, var)
                system = build_system(B, d, n)
                n_zeros = count_feasible_roots(system)
                add_row = [d n s n_zeros]
                n_eq_mat[it,:] = add_row
                it += 1
            end
        end
    end
    n_eq_mat
end

div_vec = [3 4 5 6]
hois_vec = [4 5 6]
n_sim = 3000
var = 1
data = main(div_vec, hois_vec, n_sim, var)
output_name = "dim_3_6_div_3_6_s_3000"
writedlm("../data/expected_n_roots"*output_name*".csv", data)