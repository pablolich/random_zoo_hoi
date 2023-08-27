using Random #to sample random tensors
using HomotopyContinuation #to solve systems of polynomials numerically
using LinearAlgebra #to take matrix products
using DelimitedFiles #to load and save files

"""
    buildall()

given n, d, and nsim, construct nsim kss polynomial systems of n
polynomials of degree d
"""
function buildall(n::Int32, d::Int32, nsim::Int64, vars::AbstractVector)
    allsystems = []
    #sample monomial structure
    M = monomials(vars, d)
    for i in 1:nsim

    end
end