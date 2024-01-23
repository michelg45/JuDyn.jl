"""
    diagonal(x)

Returns a diagonal matrix with the entries of the vector `x` on its diagonal. 
"""
diagonal(x) = diagm(0 => x)