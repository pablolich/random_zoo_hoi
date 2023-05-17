
function randomtensor(d, n, var, dist)
    """
    Sample entries from a tensor with dimmesions stored
    in dims from a gaussian distribution with mean 0 and sd 1
    """
    #sample tensor of size (n+1) of d dimensions
    rng = MersenneTwister()
    if dist == "normal"
        return var*randn(rng, Float64, repeat([n + 1], d)...)
    else dist == "uniform"
        return rand(rng, Float64, repeat([n + 1], d)...).-0.5
    end
end

function equation_i(T_i, vec)
    """
    Build GLV equation for species i with hois of dimension d  
    """
    #get dimension of T
    dim = ndims(T_i)
    n = size(T_i,1)
    #initialize total
    tot = 0
    if dim > 2
        #hois are more than three-way, build polynomial recursively
        for j in 1:n
            slices = repeat([Colon()], dim-1)
            tot += vec[j]*equation_i(T_i[slices...,j], vec)
        end
    elseif dim == 2
        #hois are three-way, compute quadratic form
        return dot(vec, T_i*vec)
    else
        #no hois, compute linear form
        return T_i'*vec
    end
    return tot
end

function getsystem(vars, B, d, n)
    """
    Build system of polynomials
    """
    #add a constant species
    vars = [vars; 1]
    #initialize system of ODEs as empty list
    equations = []
    #construct set of dynamic equations
    for i in 1:n
        #get tensor of interactions for species i
        B_i = B[repeat([Colon()],d-1)...,i]
        eqn = equation_i(B_i, vars)
        append!(equations, eqn)
    end
    System(equations)
end