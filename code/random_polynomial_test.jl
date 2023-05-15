#This script tests the prediction of expected number of equilibria for a 
#GLV with three-way higher order interactions

using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

function randomtensor(d, n, var, dist)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    rng = MersenneTwister()
    if dist == "normal"
        return var*randn(rng, Float64, repeat([n + 1], d)...)
    else dist == "uniform"
        return rand(rng, Float64, repeat([n + 1], d)...)
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
        #hoisa re three-way, compute quadratic form
        return dot(vec, T_i*vec)
    else
        #no hois, compute linear form
        return T_i'*vec
    end
    return tot
end

function getsystem(vars, B, d, n)
    """
    Build system of polynomials
    """
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

function randomsystem(vars, d, n)
    #initialize place holder for polynomial system
    equations = []
    #construct system
    for i in 1:n
        append!(equations, rand_poly_dist(Float64, vars, d, "normal"))
    end
    System(equations)
end

function countreal(res)
    length(real_solutions(res))
end

function filterpositive(res)
    res_pos = filter(s -> all(s .> 0), real_solutions(res))
    return res_pos
end

function main(div_vec, hois_vec, n_sim, var, dist)
    """
    Get number of positive roots for n_sim simulations
    of n species, with n running from 3 to n_max
    """
    n_hoi = length(hois_vec)
    n_div = length(div_vec)
    #preallocate matrix to store number of zeros
    global n_eq_mat = Array{Float64}(undef, 6)
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
                B = randomtensor(d, n, var, dist)
                system = getsystem(x, B, d, n)
                #show solving progress only if slow
                if n*d > 15 show = true else show = false end
                #solve system and get real solutions
                real_sols = real_solutions(solve(system, show_progress = show))
                n_real = length(real_sols)
                #get positive solutions
                pos_sols = filter(s -> all(s .> 0), real_sols)
                n_pos = length(pos_sols)
                #get largest eigenvalue of jacobian evaluated at all positive equilibria
                lambda_max_vec = zeros(Float64, n_pos)
                for i in 1:n_pos
                    j_eval_i = jacobian(system, pos_sols[i])
                    lambda_max_vec[i] = maximum(real(eigvals(j_eval_i)))
                end
                #store results
                add_rows = hcat(repeat([d n s n_real n_pos], n_pos), lambda_max_vec)
                global n_eq_mat = hcat(n_eq_mat, add_rows')
            end
        end
    end
    return n_eq_mat'
end

hoi_vec = [2 3 4 5 6]
div_vec = [3 4 5]
n_sim = 1000
var = 1
dist = "normal"
data = main(div_vec, hoi_vec, n_sim, var, dist)
#save data
max_hoi = string(maximum(hoi_vec))
max_div = string(maximum(div_vec))
output_name = "dim_"*max_hoi*"_div_"*max_hoi*"_s_"*string(n_sim)*"_"*dist
writedlm("../data/expected_n_roots_"*output_name*".csv", data)