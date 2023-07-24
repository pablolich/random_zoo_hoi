# using Optimization, ForwardDiff, Zygote, OptimizationOptimJL #optimization packages
using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically

function evaluateglv(x::Array, A::Array, B::Array, r::Array, h::Array, e::Array)
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
                append!(response, A[i,j]*x[j]/(h[j]+x[j]) - 
                B[j,i]*x[j]/(e[i] + x[i]))
            end
        end
        append!(res, r[i] - A[i,i]*x[i] + sum(response))
    end
    return res
end

function buildglvpoly(x::Array, A::Array, B::Array, r::Array, h::Array, e::Array)
    """
    function to build polynomial form of the model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        terms = []
        factors = h .+ x
        factors[i] = 1
        for j in 1:nspp
            if j != i
                append!(terms, A[i,j]*x[j]*(e[i]+x[i])*prod(factors)/(h[j]+x[j])-
                prod(factors)*B[j,i]*x[j])
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
        append!(res, (e[i]+x[i])*prod(factors)*(r[i] - A[i,i]*x[i]) + sum(terms))
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

function addnegatives(x, variables)
    """
    function to add up the negative components of solutions of system
    """
    nspp = length(variables)
    #generate parameters
    A = reshape(x, (nspp, nspp))
    H = ones(nspp, nspp)
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
        neg_comp = abs(sum(sol_i[sol_i.<0]))
        cost += neg_comp
        println(cost)
    end
    println("cost ", cost)
    return cost
end

n = 4
x0 = ones(n^2)
@var x[1:n]
f = OptimizationFunction(addnegatives, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, x0, x)
sol = Optimization.solve(prob, SimulatedAnnealing())
A = reshape(sol.u, (n, n))
H = ones(n, n)
#build system
syst = System(buildglvpoly(x, A, H))
#solve system and get real solutions
real_sols = real_solutions(HomotopyContinuation.solve(syst, show_progress=true))
#check if equilibria satisfy original model
cert_vec = certifysolutions(evaluateglv, real_sols, A, H)
cert_real_sols = real_sols[cert_vec, :]
println(cert_real_sols)