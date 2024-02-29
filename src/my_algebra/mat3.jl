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

function Mat3(n1::Vec3,n2::Vec3,n3::Vec3)

    mat = [n1.v n2.v n3.v]

    return Mat3(mat)
end