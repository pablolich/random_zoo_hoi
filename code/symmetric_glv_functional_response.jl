using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files
using Symbolics #perform analytical calculations
using Roots
using Plots

function evaluateglv(x::Array, A, r::Array, s)
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
                append!(response, A[i,j]*x[j]/(1+s*x[j]))
            end
        end
        push!(res, r[i] - x[i] - sum(response))
    end
    return res
end

function buildglvpoly(x::Array, A, r::Array, structure::BitMatrix, s)
    """
    function to build polynomial form of the model
    """
    nspp = length(x)
    res = []
    for i in 1:nspp
        terms = []
        xinteracting = x[findall(!iszero,x.*structure[i,:])]
        factors = 1 .+ s*xinteracting
        for j in 1:nspp
            if j != i
                append!(terms, A[i,j]*x[j]*prod(factors)/(1+s*x[j]))
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
        append!(res, prod(factors)*(r[i] - x[i]) - sum(terms))
    end
    return res
end

function certifysolutions(model, solutions, A, r, s)
    """
    Check if solutions are solutions to model or not
    """
    nsol = size(solutions, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(solutions[i], A, r, s)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        end
    end
    return cert_vec
end

"""
    interactionstructure(n, n_neighbor)

which species does species one interact with
"""
function interactionstructure(n, n_neighbor)
    first_row = []
    if n_neighbor == 0
        return first_row
    else
        #create vector with indices
        indices = collect(2:n)
        #reverse it
        indices_rev = reverse(indices)
        i = 1
        while i <= n_neighbor
            push!(first_row, indices[i])
            if length(first_row) == n_neighbor
                return first_row
            elseif length(first_row) > n_neighbor
                #how much longer is it
                longer = length(first_row) - n_neighbor
                for j in 0:(longer-1)
                    popat!(first_row, length(first_row)-j) 
                end
                return first_row
            else
                push!(first_row, indices_rev[i])
                i+=1
            end
        end
    end
    return first_row
end



function buildfirstrow(n, n_neighbor, interacton_strength)
    #get non-zero indices
    ind_row = interactionstructure(n, n_neighbor)
    first_row = zeros(Float64, n)
    #set those positions to desired strength
    for i in ind_row
        first_row[i] = interacton_strength
    end
    return first_row
end

function buildcirculant(first_row)
    C = []
    n = length(first_row)
    for i in 1:n
        if i == 1
            C = reshape(first_row, (1,n))
        else
            C = vcat(C, reshape(circshift(first_row, i-1), (1,n)))
        end
    end
    return C #transpose
end

function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end

function getrepeatedeigs(eigenvalues)
    #evaluate eigenvalues numerically
    numeric_eigenvals = Vector{Float64}([])
    n_eigenvals = length(eigenvalues)
    for i in 1:n_eigenvals
        push!(numeric_eigenvals, Float64(real(substitute(eigenvalues[i], Dict([a=>pi]))).val))
    end
    rounded_eigs = round.(numeric_eigenvals, digits = 3)
    rep_eigenvalues = nonunique(rounded_eigs)
    indices = findall(x->in(x, rep_eigenvalues), unique(rounded_eigs))
    return eigenvalues[indices]
end

function findcritical(lambda, parameter)
    #get real par
    real_lambda = real(lambda)
    #build function
    f_expr = build_function(real_lambda, parameter, expression=Val{false})
    f = eval(f_expr)
    critical_par = nothing
    try
        critical_par = find_zero(f, 1)
        return critical_par
    catch 
        return critical_par
    end
end 

function criticalparameters(eigenvalues, parameter)
    criticalpars = []
    n = length(eigenvalues)
    for i in 1:n
        eig = eigenvalues[i]
        push!(criticalpars, findcritical(eig, parameter))
    end
    return criticalpars[findall(!isnothing, criticalpars)]
end

#functions dealing with analytical calculations

function symmetricequilibrium(n_neighbor)
    return 1/2*(-n_neighbor*a + sqrt((n_neighbor*a)^2 + 4))
end

function jacobianelement(n_neighbor, bifurcation_par)
    return -4*bifurcation_par/(2-bifurcation_par*n_neighbor + sqrt(4+(bifurcation_par*n_neighbor)^2))^2
end

function rootsofunity(n)
    return exp(2*pi*im/n)
end

function jaceigenvals(c, n, n_neighbor)
    eigs = []
    sumind = interactionstructure(n, n_neighbor)
    for j in 0:(n-1)
        roots = 0
        for i in sumind
            roots += (rootsofunity(n))^(j*(i-1))
        end
        push!(eigs, -1 + c*roots)
    end
    return eigs
end

"""
    parametrizecritical#=  =#()

sample matrix of interactions at the bifurcation
"""
function parametrizecritical(i, j, crit_par)
    #generate first row of interaction matrix
    first_row = buildfirstrow(i, j, crit_par)
    #create matrix of interactions
    Ac = buildcirculant(first_row)
    return Ac
end

"""
    countsols(Syst)

Count number of real and positive solutions
"""
function countsols(x, A, r, s)
    interacting = A .!= 0
    syst = System(buildglvpoly(x, A, r, interacting, s))
    #solve system and get real solutions
    real_sols = real_solutions(HomotopyContinuation.solve(syst, show_progress=true))
    #check if equilibria satisfy original model
    cert_vec = certifysolutions(evaluateglv, real_sols, A, r, s)
    cert_real_sols = real_sols[cert_vec, :]
    #get number of real and positive solutions
    nsol = length(cert_real_sols)
    pos_sols = filter(s -> all(s .> 0), cert_real_sols)
    npos = length(pos_sols)
    return [nsol npos]
end

"""
    recordcriteigenvals(nmax, parameter)

function to generate data to plot all eigenvalues at critical bifurcations
"""
function recordcriteigenvals(n, parameter)
    for i in 1:n
        for j in 0:(i-1)
            println("Diversity ",i," Interactions ", j)
            #get jacobian element
            c = jacobianelement(j, parameter)
            #get eigenvalues
            eig_vec = jaceigenvals(c, i, j)
            #get only repeated ones
            eig_vec_rep = getrepeatedeigs(eig_vec)
            #get critical parameters
            par_crit = criticalparameters(eig_vec_rep, parameter)
            #how many critical points are there
            n_crit = length(par_crit)
            #initialize
            add_rows = []
            for par in 1:n_crit
                par_crit_i = par_crit[par]
                #evaluate eigenvalues at bifurcation
                numegienvals = [substitute(eig_vec[k], Dict([parameter=>par_crit_i])) for k in 1:i]
                real_part = real(numegienvals)
                imaginary_part = imag(numegienvals)
                size_connect_numb = repeat([i j par], i)
                if par == 1
                    add_rows = [size_connect_numb real_part imaginary_part]
                else
                    add_rows = vcat(add_rows, [size_connect_numb real_part imaginary_part])
                end
                #save sweep
                open("../data/eigenvalues.csv", "a") do io
                writedlm(io, add_rows, ' ')
                end
            end
        end
    end
end
#write main function

function one_sweep(n, parameter)
    for i in 1:n
        #set growthrates and variables of system
        r = ones(i)
        @var x[1:i]
        for j in (i-1):(i-1)
            println("Diversity ",i," Interactions ", j)
            #get jacobian element
            c = jacobianelement(j, parameter)
            #get eigenvalues
            eig_vec = jaceigenvals(c, i, j)
            #get only repeated ones
            eig_vec_rep = getrepeatedeigs(eig_vec)
            #get critical parameters
            par_crit = criticalparameters(eig_vec_rep, parameter)
            #how many critical points are there
            n_crit = length(par_crit)
            if n_crit == 0
                #run a randomized version of model
                first_row = buildfirstrow(i, j, rand())
                A = buildcirculant(first_row)
                add_rows = [i j countsols(x, A, r, 1) 1;
                            i j countsols(x, A, r, 0) 0]
            else
                #run model close to the each critical parametrizations
                for par in 1:n_crit
                    par_crit_i = par_crit[par]
                    #build matrix of interactions at bifurcation
                    Ac = parametrizecritical(i, j, par_crit_i)
                    interacting = Ac .!= 0
                    perturbation = rand()
                    #add small perturbation
                    A = Ac .+ 0.1*perturbation.* interacting
                    #count solutions of system
                    nsolsvec = countsols(x, A, r, 1)
                    #count solutions of control syste
                    nsolsvec_control = countsols(x, A, r, 0)
                    if par == 1
                        #add treatement
                        add_rows = [i j nsolsvec 1;
                                    i j nsolsvec_control 0]
                        #add control
                    else
                        add_rows = vcat(add_rows, [i j nsolsvec 1;
                                                   i j nsolsvec_control 0])
                    end
                end
            end
            #save sweep
            open("../data/parameter_sweeps_glvfunc.csv", "a") do io
                writedlm(io, add_rows, ' ')
            end
        end
    end
end

function sweep_many(n_sweeps, max_size, parameter)
    count = 0
    while count < n_sweeps
        one_sweep(max_size, parameter)
        count+=1
    end
end

@variables a
#recordcriteigenvals(15, a)
sweep_many(1000, 7, a)