using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

function evaluateglv(x::Array, A::Array, H::Array)
    """
    function to evaluate original model
    """
    nspp = length(x)
    res = ones(nspp)
    for i in 1:nspp
        response = zeros(nspp)
        for j in 1:nspp
            if j == i
                response[j] = 0
            else 
                response[j] = A[i,j]*x[j]/(1+H[i,j]*x[j])
            end
        end
        res[i] = 1 - A[i,i]*x[i] - sum(response)
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
        terms = zeros(nspp)
        Hrowi = H[i,:]
        factors = 1 .+ Hrowi.*x
        factors[i] = 1
        for j in 1:nspp
            if j != i
                factors[j] = 1
                terms[j] = A[i,j]*x[j]*prod(factors)
            end
        append!(res, prod(factors)*(1 - A[i,i]*x[i]) - sum(terms))
        end
    end
    return res
end

function certifysolutions(model, solutions, x, A, H)
    nsol = size(solutions, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(x, A, H)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        else push!(cert_vec, i)
        end
    end
    return cert_vec
end

function onesweep(nspp_vec)
    for i in nspp_vec
        #generate parameters
        A = rand(Float64, (i, i))
        H = rand(Float64, (i, i))
        @var x[1:i]
        #build system
        syst = System(buildglvpoly(x, A, H))
        #solve system and get real solutions
        real_sols = real_solutions(solve(syst, show_progress=true))
        #check if equilibria satisfy original model
        cert_vec = certifysolutions(evaluateglv, real_sols, x, A, H)
        cert_real_sols = real_sols[cert_vec, :]
        #get number of real and positive solutions
        nsol = length(cert_real_sols)
        pos_sols = filter(s -> all(s .> 0), cert_real_sols)
        npos = length(pos_sols)
        add_rows = [i nsol npos]
        #save sweep
        open("../data/parameter_sweeps_glvfunc.csv", "a") do io
            writedlm(io, add_rows, ' ')
        end    
    end
end

function manysweeps(nspp_vec, nsim)
    s=0
    while s < nsim
        onesweep(nspp_vec)
        s+=1
    end
end

manysweeps([1,2,3,4,5], 2)

#Plant an equilibrium

#Come up with parameters that minimize negative solutions and see if I maintain number of 
#real