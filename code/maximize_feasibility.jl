using Optimization, ForwardDiff, Zygote, OptimizationOptimJL #optimization packages
using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically

function evaluateglv(x::Array, A::Array, H::Array)
    """
    function to evaluate original model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        response = []
        for j in 1:nspp
            if j == i
                append!(response, 0)
            else 
                append!(response, A[i,j]*x[j]/(1+H[i,j]*x[j]))
            end
        end
        append!(res, 1 - A[i,i]*x[i] - sum(response))
    end
    return res
end

function buildglvpoly(x::Array, A::Array, H::Array)
    """
    function to build polynomial form of the model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        terms = []
        Hrowi = H[i,:]
        factors = 1 .+ Hrowi.*x
        factors[i] = 1
        for j in 1:nspp
            if j != i
                append!(terms, A[i,j]*x[j]*prod(factors)/(1+H[i,j]*x[j]))
            end
            #try/catch to handle the particularly annoying case of one species
            try
                empty = !any(terms)
                if empty
                    append!(terms, 0)
                end
            catch
            end
        end
        append!(res, prod(factors)*(1 - A[i,i]*x[i]) - sum(terms))
    end
    return res
end

function certifysolutions(model, solutions, A, H)
    nsol = size(solutions, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(solutions[i], A, H)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        end
    end
    return cert_vec
end

mutable struct Parameters
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}
    growthalleles::Vector{Float64}
    growthrates::Vector{Float64}
end

function addnegatives(x, variables)
    """
    function to add up the negative components of solutions of system
    """
    nspp = length(variables)
    #generate parameters
    A = reshape(x, (nspp, nspp))
    H = ones(i,i)
    #build system
    syst = System(buildglvpoly(variables, A, H))
    #solve system and get real solutions
    real_sols = real_solutions(solve(syst, show_progress=true))
    #check if equilibria satisfy original model
    cert_vec = certifysolutions(evaluateglv, real_sols, A, H)
    cert_real_sols = real_sols[cert_vec, :]
    #get number of real solutions
    nsol = length(cert_real_sols)
    cost = 0
    for sol in 1:nsol
        sol_i = cert_real_sols[sol]
        neg_comp = abs(sum(cert_real_sols[cert_real_sols.<0]))
        cost += neg_comp
    end
    return cost
end

n = 2
x0 = ones(n^2)
@var x[1:n]
f = OptimizationFunction(addnegatives, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, x0, x)
sol = solve(prob, NelderMead())