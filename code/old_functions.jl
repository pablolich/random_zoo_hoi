
function multinomialcoeff(n, kvec)
    """
    Compute multinomial coefficient
    """
    num = factorial(big(n))
    den = prod([factorial(big(i)) for i in kvec])
    return num/den
end

function rand_poly_dist(T, 
    vars::AbstractVector, 
    d::Integer,
    dist::String; 
    homogeneous::Bool = false
)
    """
    Create a random dense polynomial of degree `d` in the given variables `variables`.
    Each coefficient is sampled independently from a Normal(0, var) or a Uniform(-0.5,0.5).
    """
    M = monomials(vars, d; affine = !homogeneous)
    n_terms = length(M)
    coefficient_list = zeros(Float64, n_terms)
    if dist == "normal"
        for i in 1:n_terms
            monomial_i = M[i]
            #get exponents of each variable in monomial i
            exponents, coeffs = exponents_coefficients(monomial_i, vars)
            monomial_degree = degree(monomial_i)
            #compute the variance of ith coefficient using mulitinomial coefficient
            variance = multinomialcoeff(monomial_degree, exponents)
            #sample ith coefficient from a gaussian with computed variance
            coefficient_list[i] = variance*randn(T)
        end
        sum(coefficient_list .* M)
    else dist == "uniform"
        sum((rand(T, length(M)).-0.5).*M)
    end
end

function randomsystem(vars, d, n, dist)
    """
    Build system of random polynomials
    """
    #initialize place holder for polynomial system
    equations = []
    #construct system
    for i in 1:n
        append!(equations, rand_poly_dist(Float64, vars, d, dist))
    end
    System(equations)
end
