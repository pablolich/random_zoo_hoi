using SciMLBase
using LinearAlgebra
using HomotopyContinuation
using Random
using DelimitedFiles
using DifferentialEquations
using Plots
using Statistics


"""
    glvfunc!(dx,x,p,t)

lotka-volterra with type II functional response
"""
function glvfunc!(dx,x,p,t)
    #parse parameters
    r, A = p
    dx[1] = x[1]*(r[1]-A[1,1]*x[1] - A[1,2]*x[2]/(1+x[2]))
    dx[2] = x[2]*(r[2] - A[2,2]*x[2] - A[2,1]*x[1]/(1+x[1]))
    return dx
end

"""
    glvfuncpoly!(dx,x,p,t)

lotka-volterra with hois coming from non-linearities
"""
function glvfuncpoly!(dx,x,p,t)
    r, A = p
    dx[1] = x[1]*((r[1]-A[1,1]*x[1])*(1+x[2]) - A[1,2]*x[2])
    dx[2] = x[2]*((r[2] - A[2,2]*x[2])*(1+x[1]) - A[2,1]*x[1])
end

"""
    glvext!(dy, y, p, t)

a glv model representing a non-linear system via glv-trick
"""
function glvext!(dy, y, p, t)
    s, W, T, n = p
    #get variables of original system
    x = last(y, n)
    #ensure variable constraints
    y = exp.(transpose(T)*log.(abs.(x))) #take absoltue values of x ot avoid errors
    #run extended glv equations
    #Dy = Diagonal(y)
    for i in 1:length(y)
        dy[i] =  y[i] * (s[i] + dot(W[i,:], y))    
    end
    return dy
end

"""
    glvhoi!(dx, x, p, t)

glv-hoi model for integration
"""
function glvhoi!(dx, x, p, t)
    r, A, B = p
    n = length(r)
    for i in 1:n
        pairs = 0
        triplets = 0
        for j in 1:n
            pairs += A[i,j]*x[j]
            for k in 1:n
                triplets += B[i,j,k]*x[j]*x[k]
            end
        end
        dx[i] = x[i]*(r[i]+ pairs + triplets)
    end
    return dx
end

"""
    storerow(i::Int64, storemat, tostore::Matrix{Int64})

storing results in a matrix row wise
"""
function storerow(i::Int64, storemat, tostore)
    if i == 1
        storemat = tostore
    else
        storemat = vcat(storemat, tostore)
    end
    return storemat
end
"""
    storerow(i::Int64, storemat, tostore::Matrix{Int64})

storing results in a matrix column wise
"""
function storecol(i::Int64, storemat, tostore)
    if i == 1
        storemat = tostore
    else
        storemat = hcat(storemat, tostore)
    end
    return storemat
end

"""
    variance(d::Int64, j_vec::Array)

calculate variance of coefficitent whose associated monomial has powers j_vec
"""
function variance(d::Int64, j_vec::Array)
    j_vec_expanded = [j_vec; (d - sum(j_vec))]
    multinomial(d, j_vec_expanded)
end

"""
    multinomialcoeff(n, kvec)

compute multinomial coefficient
"""
function multinomial(n::Int64, kvec::Array)
    num = factorial(n)
    den = prod([factorial(i) for i in kvec])
    num/den
end

"""
    buildpoly(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector)

given all the monomials, sample random coefficients for each and sum them to build a polynomial
"""
function buildpoly(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, d::Int64, rng::AbstractRNG)
    #initialize list of polynomial coefficients
    coefficient_list = zeros(Float16, nmon)
    for i in 1:nmon
        monomial_i = allmonomials[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, vars; expanded = true) #expanded reduces time
        #compute the variance of ith coefficient using mulitinomial coefficient
        var = variance(d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(var)*randn(rng, Float16)
    end
    sum(coefficient_list .* allmonomials)
end

"""
    buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64, rng::AbstractRNG)

construct a kss system of n polynomials of degree d
"""
function buildsystem(allmonomials::Vector{Expression}, nmon::Int64, vars::AbstractVector, 
    n::Int64, d::Int64, rng::AbstractRNG)
    equations = []
    #loop over number of equations in each system
    for j in 1:n
        append!(equations, buildpoly(allmonomials, nmon, vars, d, rng))
    end
    System(equations)
end

"""
    coeffstosyst(x, coeffs)

transform glv-hoi coefficients into multivariate polynomial system.
"""
function coeffstosyst(x, coeffs)
    r, A, B = coeffs
    n = length(r)
    dy = []
    for i in 1:n
        pairs = 0
        triplets = 0
        for j in 1:n
            pairs += A[i,j]*x[j]
            for k in 1:n
                triplets += B[i,j,k]*x[j]*x[k]
            end
        end
        row_i = r[i]+ pairs + triplets
        dy = storerow(i, dy, row_i)
    end
    return System(dy)
end

"""
    isoscillating(solution::ODESolution, n::Int64)

check if average of oscillations converges to constant solution
"""
function  isoscillating(solution::ODESolution, n::Int64)
    #get dense solution at 10 different equispaced time intervals
    t_dense = collect(range(solution.t[1], solution.t[end], 1000))
    t_intervals = collect(range(solution.t[1], solution.t[end], 10))
    soldense = solution(t_dense)
    #get abundances over the last tenth corresponding only to the real variables
    indslast = findall(t_dense .> t_intervals[9])
    solutiontail = soldense[end-n+1:end, indslast]
    npoints = length(indslast)
    #calculate moving average
    nvar, nt = size(solutiontail)
    moving_av = []
    for i in 1:nt
        #get a solution window
        sol_window = solutiontail[:,1:i]
        #calculae average for such window
        av_sol = mean(sol_window, dims=2)
        moving_av = storecol(i, moving_av, av_sol)
    end
    #check if the average becomes constant
    differences = diff(moving_av, dims = 2)
    #get the sum of the absoltue values of those differences for each time series
    total_differences = sum(abs.(differences), dims = 2)
    #impose that the total difference for all species is small
    if all(total_differences .< 1e-4*npoints)
        return true
    else
        return false
    end
end

"""
    isconstant(solution::ODESolution, n::Int64)

check if solution is constant
"""
function isconstant(solution::ODESolution, n::Int64)
    #get dense vector of equidistant timepoints
    t_dense = collect(range(solution.t[1], solution.t[end], 1000))
    #get sparse vector of equidistant timepoints
    t_intervals = collect(range(solution.t[1], solution.t[end], 10))
    #evaluate solution at dense vector of time points 
    soldense = solution(t_dense)
    #get abundances of the last tenth of the dense interval corresponding only to real species
    indslast = findall(t_dense .> t_intervals[9])
    npoints = length(indslast)
    abundanceslast = soldense[end-n+1:end, indslast]
    #get differences between consecutive abundances of each species
    differences = diff(abundanceslast, dims = 2)
    #get the sum of the absoltue values of those differences for each time series
    total_differences = sum(abs.(differences), dims = 2)
    #impose that the total difference for all species is small
    if all(total_differences .< 1e-4*npoints)
        return true
    else
        return false
    end
end

"""
    ispersistent(solution::ODESolution, n::Int64)

check if solution is persistent (steady state or bounded oscillations)
"""
function ispersistent(solution::ODESolution, n::Int64)
        if isoscillating(solution, n)
            println("solution is oscillating")
            return true
        elseif isconstant(solution, n)
            println("soluttion is constant")
            return true
        else
            println("solution did not converge")
            return false
        end
end

"""
    isdivergent(solution::ODESolution, maxtol::Float64, n::Int64)

check if solution is diverging
"""
function isdivergent(solution::ODESolution, maxtol::Float64, n::Int64)
    #get only abundances of the real species
    endstate = last(solution[end], n)
    #check if the end state is divergent
    if all(abs.(endstate) .< maxtol)
        println("solution did not diverge")
        return true    
    else 
        println("solution diverged")
        return false
    end
end

"""
    getstablediversity(func, initial::Vector{Float64}, tspan::Tuple, parameters::Tuple, n::Int64)

integrate glv dynamics given initial conditions, parameters, and time span, until persistent state is reached
"""
function getstablediversity(func, initial::Vector{Float64}, tspan::Tuple, parameters::Tuple, n::Int64)
    #create ODE problem
    problem = ODEProblem(func, initial, tspan, parameters)
    #solve it
    sol = DifferentialEquations.solve(problem, Tsit5())
    if isdivergent(sol, 1e6, n)
        return false
    else
        #check if a persistent state has been reached
        persistent = ispersistent(sol, n)
        while !persistent
            println("re-integrating...")
            #set new initial state as last state of the previous integration
            initial = sol[end]
            #zero out extinct species
            initial[findall(initial.>1e-6)] .= 0
            #integrate again
            problem = ODEProblem(func, initial, tspan, parameters)
            sol = DifferentialEquations.solve(problem, Tsit5())
            #check if solution diverges in subsequent integrations
            if isdivergent(sol, 1e6, n)
                return false
            else
                persistent = ispersistent(sol, n)
            end
        end
    end
    endstate = sol[end]
    diversity = length(findall(endstate.>1e-6))
    return diversity
end

"""
    getexpcoeffsyst(system::System, n::Int64, d::Int64, vars::AbstractVector)

Get matrices of exponents and coefficients of the system
"""
function getexpcoeffsyst(system::System, n::Int64, d::Int64, vars::AbstractVector)
    #initialize coeffficients matrix
    O = []
    r = []
    for i in 1:n
        #get equation i's exponents and coefficients
        exps, coeffs = exponents_coefficients(system.expressions[i], vars)
        #leave 0 esxponents out
        coeffs0 = coeffs[1:(end-1)]
        coeffsr = coeffs[end]
        #store
        O = storerow(i, O, reshape(coeffs0, (1, length(coeffs0))))
        r = storerow(i, r, coeffsr)
    end
    exps, coeffs = exponents_coefficients(system.expressions[1], vars)
    T = exps[:, 1:(end-1)]
    return T, O, r
end

"""
    integratekss(n::Int64, d::Int64, vars::AbstractVector, rng::AbstractRNG)

build and integrate a kss system of n polynomials of degree d
"""
function integratekss(n::Int64, d::Int64, vars::AbstractVector, rng::AbstractRNG)
    #build system
    allmon = monomials(vars, d)
    nmon = length(allmon)
    syst = buildsystem(allmon, nmon, vars, n, d, rng)
    #get parameters and initial condition of corresponding glv system
    T, O, r = getexpcoeffsyst(syst, n, d, vars)
    W = transpose(T)*O #extended matrix of interactions
    s = transpose(T)*r #extended growth rates
    pglv = (s, W, T, n)
    initial = ones(size(T, 1))
    #integrate dynamics
    div = getstablediversity(glvext!, initial, (0, 1e6), pglv, n)
    return div
end

"""
    sweepntimes(max_n::Int64, max_d::Int64, nsweeps::Int64)

generate data§ from integrations
"""
function sweepntimes(max_n::Int64, max_d::Int64, nsweeps::Int64, seed::Int64)
    parameters = [(x, y) for x in 1:max_n, y in 1:max_d]
    n_pairs = length(parameters)
    #initialize random generator
    rng = MersenneTwister(seed)
    #initialize storing 
    results = []
    for n_d in 1:n_pairs
        n = parameters[n_d][1]
        d = parameters[n_d][2]
        @var x[1:n]
        div = false
        for sim in 1:nsweeps
            println("System size ", n, " System degree ", d, " Simulation ", sim)
            while !div
                div = integratekss(n, d, x, rng)
            end
            append!(results,[n d div])
        end
    end
    open("../data/kss_simulations.csv", "a") do io
        writedlm(io, results', ' ')
    end
end

sweepntimes(3, 3, 2, 1)

"""
    integrateitmescale()

generate data for timescale plots
"""
function integratetimescale()
    x0 = [0.5;1.5]
    r = [1.0;1.0]
    A = [1.0 0.1;0.2 1.0]
    p = (r, A)
    tspan = (0.0, 10.0)
    prob = ODEProblem(glvfunc!, x0, tspan, p)
    sol = DifferentialEquations.solve(prob, Tsit5())
    prob = ODEProblem(glvfuncpoly!, x0, tspan, p)
    solpoly = DifferentialEquations.solve(prob, Tsit5())
    #get data frame with dense solutions
    t_dense = collect(range(0, stop=10, length=1000))
    solt = sol(t_dense)
    soltpoly = solpoly(t_dense)
    solsall = [t_dense solt[1,:] solt[2,:] soltpoly[1,:] soltpoly[2,:]]
    open("../data/timescaleseparation.csv", "a") do io
        writedlm(io, solsall, ' ')
    end
end

"""
    testglvtrick()

test that the glv trick gives me the same results as integrating the dynamics directly
"""
function testglvtrick()
    #set variables and integration time window
    @var z[1:3]
    tspan = (0.0, 20.0)
    #initial conditions and parameters of glv with hoi
    x0 = [0.5;1.5;2.0]
    r = [1.0;1.0;1.0]
    A = -[1.0 0.1 0.1;0.2 1.0 0.2; 0.3 0.3 1.0]
    B = -rand(3,3,3)
    phoi = (r, A, B)
    #integrate normal glv hoi
    prob = ODEProblem(glvhoi!, x0, tspan, phoi)
    solhoi = DifferentialEquations.solve(prob, Tsit5())
    stationaryhoi = solhoi[end]

    #transform hoi parameters to polynomial system
    syst = coeffstosyst(z, phoi)
    #get parameters and initial conditions for equivalent higher dimensional glv
    T, O, r = getexpcoeffsyst(syst, 3, 3, z)
    W = transpose(T)*O #matrix of interactions
    s = transpose(T)*r #growth rates
    pglv = (s, W, T, 3)
    y0 = exp.(transpose(T)*log.(x0))
    #integrate extended glv
    prob = ODEProblem(glvext!, y0, tspan, pglv)
    sol = DifferentialEquations.solve(prob, Tsit5())
    stationarytrick = sol[7:9,end]
    println("Solution from brute force integration: ", stationaryhoi)
    println("Solution from glv-trick: ", stationarytrick)
end