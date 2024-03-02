"""
    diagonal

Returns a diagonal matrix with the entries of the vector `x::Vector{Float64}` on its diagonal. 

Calling sequence:

    d = diagonal(x)

"""
diagonal(x::Vector{Float64}) = diagm(0 => x)