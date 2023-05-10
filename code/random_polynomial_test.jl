#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles #to load and save files

function randomtensor(d, n, var)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    rng = MersenneTwister()
    #mean 0, variance var
    var*randn(rng, Float64, repeat([n + 1], d)...)
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
        #hoisa re three-way, compute quadratic form
        return dot(vec, T_i*vec)
    else
        #no hois, compute linear form
        return T_i'*vec
    end
    return tot
end

function getsystem(B, d, n)
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
        #get tensor of interactions for species i
        B_i = B[repeat([Colon()],d-1)...,i]
        eqn = equation_i(B_i, x)
        append!(equations, eqn)
    end
    System(equations)
end

function countpositive(system)
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
                B = randomtensor(d, n, var)
                system = getsystem(B, d, n)
                n_zeros = countpositive(system)
                add_row = [d n s n_zeros]
                n_eq_mat[it,:] = add_row
                it += 1
            end
        end
    end
    n_eq_mat
end

hoi_vec = [2 3 4 5 6]
div_vec = [3 4 5]
n_sim = 2000
var = 1
data = main(div_vec, hoi_vec, n_sim, var)
#save data
max_hoi = string(maximum(hoi_vec))
max_div = string(maximum(div_vec))
output_name = "dim_"*max_hoi*"_div_"*max_hoi*"_s_"*string(n_sim)
writedlm("../data/expected_n_roots"*output_name*".csv", data)