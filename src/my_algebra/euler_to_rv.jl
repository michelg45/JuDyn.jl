"""
euler_to_RV

    function that computes the rotation vector resulting from a set of Euler angles

        input: phi, theta, psi::Float64     Euler angles (in radians)
        output: rv::RV3                     rotation vector

        calling sequence: rv = euler_to_RV(phi,theta,psi)

        Example:

        psi = 0.0
        0.000000e+00
        theta = pi/9.0
        3.490659e-01
        phi = 0.0
        0.000000e+00
        rv = euler_to_RV(phi,theta,psi)
        RV3([3.490659e-01, 0.000000e+00, 0.000000e+00])

"""
function euler_to_RV(phi::Float64,theta::Float64,psi::Float64)


    R=rot(psi,3)*rot(theta,1)*rot(phi,3)

    return  invrot(R)

    end
