"""
    Mat3

        mutable structure hosting a 3x3 matrix (e.g. rotation operator, tangent operator)

        Examples

```@example
        phi = RV3(0.0, 2.0, -1.0)
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])

        A = tilde(phi)
        Mat3([0.000000e+00 1.000000e+00 2.000000e+00; 
        -1.000000e+00 0.000000e+00 -0.000000e+00;
        -2.000000e+00 0.000000e+00 0.000000e+00])

        R = rot(phi)
        Mat3([-6.172729e-01 3.518449e-01 7.036898e-01; 
        -3.518449e-01 6.765454e-01 -6.469092e-01; 
        -7.036898e-01 -6.469092e-01 -2.938183e-01])

        T = tang(phi)
        Mat3([3.518449e-01 -3.234546e-01 -6.469092e-01; `
        3.234546e-01 8.703690e-01 -2.592620e-01; 
        6.469092e-01 -2.592620e-01 4.814759e-01])

```

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