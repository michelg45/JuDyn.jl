

"""
    RV3(a::RV3,b::RV3)

        function performing the product of successive rotations in terms of rotation parmeters.

            calling sequence: c = RV3(a,b)

            Example

            phi = RV3(0.0, 2.0, -1.0)
            RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])
            psi = RV3(5.0, 2.0, +1.0)
            RV3([5.000000e+00, 2.000000e+00, 1.000000e+00])
            theta = RV3(phi,psi)
            RV3([-6.605006e-01, 1.995766e+00, -2.754112e-01])

"""
function RV3(a::RV3,b::RV3)


    return invrot(rot(a)*rot(b))

end


