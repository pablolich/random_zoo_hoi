using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

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
        factors = h + x
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

function certifysolutions(model, solutions, A, B, r, h, e)
    nsol = size(solutions, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(solutions[i], A, B, r, h, e)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        end
    end
    return cert_vec
end

function onesweep(nspp_vec)
    for i in nspp_vec
        #generate parameters
        A = 0.6*ones((i,i)) .+ 0.01*rand((i, i))
        B = 0.2*ones((i, i)).+ 0.005*rand((i, i))
        r = 0.3*ones(i) + 0.01*rand(i)
        h = 0.3*ones(i) + 0.01*rand(i)
        e = 0.3*ones(i) + 0.01*rand(i)
        for k in 1:i A[k,k] = 0.01 end
        #A = rand(Float64, (i, i))
        #B = rand(Float64, (i, i))
        #r = rand(Float64, i)
        #h = rand(Float64, i)
        #e = rand(Float64, i)
        @var x[1:i]
        #build system
        syst = System(buildglvpoly(x, A, B, r, h, e))
        #solve system and get real solutions
        real_sols = real_solutions(solve(syst, show_progress=true))
        #check if equilibria satisfy original model
        cert_vec = certifysolutions(evaluateglv, real_sols, A, B, r, h, e)
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

manysweeps([1, 2, 3, 4, 5, 6, 7, 8], 100)