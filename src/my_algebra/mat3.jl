"""
    Mat3

mutable structure hosting a 3x3 matrix (e.g. rotation operator, tangent operator)

"""

mutable struct Mat3
    mat::Array{Float64,2}
end

function Mat3()
    return Mat3(zeros(3,3))
end
