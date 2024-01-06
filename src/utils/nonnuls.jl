"""
    nonnuls

        function returninng a vector with the non-zero elements of a vector ::Vector(Int).

            calling sequence= b = nonnuls(a)
"""
nonnuls(a::Vector{Int}) = a[a .> 0]
