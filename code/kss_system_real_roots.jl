using HomotopyContinuation #to solve systems of polynomials numerically
using Statistics

function multinomialcoeff(n, kvec)
    """
    Compute multinomial coefficient
    """
    num = factorial(n)
    den = prod([factorial(i) for i in kvec])
    return num/(den*factorial(n-sum(kvec)))
end

function rand_poly_dist(T, 
    vars::AbstractVector, 
    d::Integer,
    homogeneous::Bool = false
)
    """
    Create a random dense polynomial of degree `d` in the given variables `variables`.
    Each coefficient is sampled independently from a Normal(0, var) or a Uniform(-0.5,0.5).
    """
    M = monomials(vars, d; affine = !homogeneous)
    n_terms = length(M)
    coefficient_list = zeros(Float64, n_terms)
    for i in 1:n_terms
        monomial_i = M[i]
        #get exponents of each variable in monomial i
        exponents, coeffs = exponents_coefficients(monomial_i, x)
        #monomial_degree = degree(monomial_i)
        #compute the variance of ith coefficient using mulitinomial coefficient
        variance = multinomialcoeff(d, exponents)
        #sample ith coefficient from a gaussian with computed variance
        coefficient_list[i] = sqrt(variance)*randn(T)
    end
    #println(sum(coefficient_list .* M))
    sum(coefficient_list .* M)
end

function randomsystem(vars, d, n)
    """
    Build system of random polynomials
    """
    
    #initialize place holder for polynomial system
    equations = []
    #construct system
    for i in 1:n
        append!(equations, rand_poly_dist(Float64, vars, d))
    end
    System(equations)
end

function main(n, d, n_sim)
    """
    Run bulk of simulations
    """
    #variables 
    @var x[1:n]
    nsol = []
    expected = sqrt((d)^(n))
    for s in 1:n_sim
        syst = randomsystem(x, d, n)
        #solve system and get real solutions
        real_sols = real_solutions(solve(syst))
        append!(nsol, length(real_sols))
        println("expected ", expected," ", mean(nsol), " observed")

    end
end
 
main(2, 3, 1000)