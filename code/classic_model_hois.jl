using Random #to sample random parameters
using Distributions #sample from special distributons
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

function sumterm(C, kappa, mu, nres, i, solution)    
    #useful vector
    Cx = C*solution
    terms = []
    for j in 1:nres
        term_j = C[j,i]*kappa[j]/(mu[j] + Cx[j])
        append!(terms, term_j)
    end
    return sum(terms)
end

function testsolution(solution, C, kappa, mu, d)
    nres, nspp = size(C)
    dxstardt = []
    for i in 1:nspp
        term1_i = sumterm(C, kappa, mu, nres, i, solution)
        dxstardt_i = term1_i - d[i]
        append!(dxstardt, dxstardt_i)
    end
    return dxstardt
end

function generatepars(distribution, nspp, nres)
    C = rand(distribution, nres, nspp)
    mu = rand(distribution, nres, 1)
    kappa = rand(distribution, nres, 1)
    death = rand(distribution, nspp, 1)
    return C, mu, kappa, death
end

function test_pars(nspp, nres)
    C = ones(nres, nspp)
    mu = 0.5*ones(nres, 1)
    kappa = 4*ones(nres, 1)
    death = 2*ones(nspp, 1)
    return C, mu, kappa, death
end

function product1out(nres, j, factors)
    #indices to keep
    keep = deleteat!(collect(1:nres), j)
    #take out jth term in the product
    prod_poly = prod(factors[keep])
    return prod_poly
end

function sum_terms_j(nres, factors, C_mat, kappa_vec, i)
    sum_j_terms = []
    for j in 1:nres
        append!(sum_j_terms, product1out(nres, j, factors)*C_mat[j,i]*kappa_vec[j])
    end
    return sum(sum_j_terms)
end

function buildequation(nres, i, factors, C_mat, kappa_vec, death_vec)
    term1 = sum_terms_j(nres, factors, C_mat, kappa_vec, i)
    term2 = death_vec[i]*prod(factors)
    return term1 - term2
end

function buildsystem(C_mat, mu_vec, kappa_vec, death_vec, variables)
    """
    Writes the polynomials S_k(t)
    C_mat [m x n]
    mu_vec [m x 1]
    kappa_vec [m x 1]
    death_vec [n x 1]
    """
    #number of species and resources
    nres, nspp = size(C_mat)
    #build all vectors that will be used
    factors = mu_vec + C_mat*variables #elements of the product
    equations = []
    for i in 1:nspp
        equation_i = buildequation(nres, i, factors, C_mat, kappa_vec, death_vec)
        append!(equations, equation_i)
    end
    return System(equations)
end

function main(iterator, nsim)
    #preallocate matrix to store results
    n_eq_mat = zeros(Float64, 1, 4)
    #loop over species-resource pairs
    for n_m in iterator
        #get number of species and resources in the pair
        n=n_m[1]
        m=n_m[2]
        #declare dynamic variables
        @var x[1:n]
        #loop over simulations
        println("Species ", n)
        println("Resources ", m)
        for s in 1:nsim
            #sample parameters
            #C, mu, kappa, d = generatepars(LogNormal(), n, m)
            C, mu, kappa, d = test_pars(n, m)
            syst = buildsystem(C, mu, kappa, d, x)
            #solve system and get real solutions
            real_sols = real_solutions(solve(syst))
            #test solution
            nsol = length(real_sols)
            #get positive solutions
            pos_sols = filter(s -> all(s .> 0), real_sols)
            npos = length(pos_sols)
            #test solutions
            if npos > 0
                for i in 1:npos
                    test = testsolution(pos_sols[i], C, kappa, mu, d)
                    println(test)
                end
            end
            #store results
            add_rows = [n m nsol npos]
            n_eq_mat = vcat(n_eq_mat, add_rows)
        end
    end
    return n_eq_mat
end

nspp_max = 3
nres_max = 3
#set parameters for simulations
spp_res_pairs = [(x, y) for x in 1:nspp_max for y in 1:nres_max if y >= x]
nsim = 100
#run simulations
data = main(spp_res_pairs, nsim)
#save results
writedlm("../data/cr_feasibility.csv", data)