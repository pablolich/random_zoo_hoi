using SciMLBase
using LinearAlgebra
using HomotopyContinuation
using Random
using DelimitedFiles
using DifferentialEquations


"""
    glv!(dy, y, p, t)

a glv model
"""
function glv!(dy, y, p, t)
    s, W = p
    Dy = Diagonal(y)
    dy =  Dy * (s - W*y)    
end

"""
    glvhoi!(dx, x, p, t)

glv-hoi model for integration
"""
function glvhoi!(dx, x, p, t)
    r, A, B = p
    n = lengh(r)
    for i in 1:n
        pairs = []
        triplets = []
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
        row_i = x[i]*(r[i]+ pairs + triplets)
        dy = storerow(i, dy, row_i)
    end
    return System(dy)
end

# set variables and integration time window
@var x[1:3]
tspan = (0.0, 10.0)
#initial conditions and parameters of glv with hoi
x0 = [0.5;1.5;2.0]
r = [1.0;1.0;1.0]
A = [1.0 0.1 0.1;0.2 1.0 0.2; 0.3 0.3 1.0]
B = rand(3,3,3)
phoi = (r, A, B)
#integrate normal glv hoi
prob = ODEProblem(glvhoi!, x0, tspan, phoi)
solhoi = DifferentialEquations.solve(prob, Tsit5())

#transform hoi parameters to polynomial system
syst = coeffstosyst(x, phoi)
#get parameters and initial conditions for equivalent higher dimensional glv
T, O = getexpcoeffsyst(syst, 3, 3, x)
W = transpose(T)*O #matrix of interactions
s = transpose(T)*r #growth rates
y0 = exp(transpose(T)*log(x0))
#integrate extended glv
prob = ODEProblem(glv!, x0, tspan, pglv)
sol = DifferentialEquations.solve(prob, Tsit5())

#get data frame with dense solutions
t_dense = collect(range(0, stop=10, length=1000))
solt = sol(t_dense)
soltpoly = solpoly(t_dense)

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
    integrateitmescale()

generate data for timescale plots
"""
function integrateitmescale()
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
    storerow(i::Int64, storemat, tostore::Matrix{Int64})

storing results in a matrix
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
    getexpcoeffsyst(system::System, n::Int64, d::Int64, vars::AbstractVector)

Get matrices of exponents and coefficients of the system
"""
function getexpcoeffsyst(system::System, n::Int64, d::Int64, vars::AbstractVector)
    #initialize coeffficients matrix
    O = []
    for i in 1:n
        #get equation i's exponents and coefficients
        exps, coeffs = exponents_coefficients(system.expressions[i], vars)
        #leave 0 exponents out
        coeffs0 = coeffs[1:(end-1)]
        coeffsr = coeffs[end]
        #store
        O = storerow(i, O, reshape(coeffs0, (1, length(coeffs0))))
    end
    T = exps[:, 1:(end-1)]
    return T, O
end

function main(vars::AbstractVector, n::Int64, d::Int64, rng::AbstractRNG)
    allmon = monomials(vars, d)
    #build system
    syst = buildsystem(allmon, length(allmon), vars, n, d, rng)
    #build matrix W
    W = getparsglv(syst, n, d, vars)
end

seed = 1
rng = MersenneTwister(seed)
n = 2
d = 2
@var x[1:n]
#y0 = 