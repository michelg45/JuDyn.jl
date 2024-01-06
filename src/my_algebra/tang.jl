"""
    tang(psi::RV3,a::Vec3)
        function computing the product of a vector a::Vec3 
        by a tangent operator with Cartesian vector psi::RV3

        calling sequence: b = tang(psi,a)

        Example

        a = Vec3(1.0,0.0,0.0)
        Vec3([1.000000e+00, 0.000000e+00, 0.000000e+00])
        psi = RV3(0.0,2.1,0.0)
        RV3([0.000000e+00, 2.100000e+00, 0.000000e+00])
        b = tang(psi,a)
        Vec3([4.110521e-01, 0.000000e+00, 7.165934e-01])



"""
function tang(psi::RV3,a::Vec3)
    psi = psi.v
    a = a.v
    PSI2 = psi'*psi
    PSI=sqrt(PSI2)
    PSI > 1.0e-2 ? f_1 = (cos(PSI)-1.0)/PSI2 : f_1 = -0.5 + PSI2/24.0 - PSI2^2/720.0
    PSI > 1.0e-4 ? f_2 = (1.0-sin(PSI)/PSI)/PSI^2 : f_2 = 1.0/6.0 - PSI2/120.0 + PSI2^2/5040.0
    b = a + f_1*cross(psi,a) + f_2*cross(psi,cross(psi,a))
    return Vec3(b)
end

"""
    tang(psi::RV3)
        function computing the tangent operator T::Mat3
        from the Cartesian rotation vector psi::RV3 

        calling sequence: T = tang(psi)

        psi = RV3(0.0,2.1,0.0)
        RV3([0.000000e+00, 2.100000e+00, 0.000000e+00])
        T = tang(psi)
        Mat3([4.110521e-01 0.000000e+00 -7.165934e-01; 
        0.000000e+00 1.000000e+00 0.000000e+00; 
        7.165934e-01 0.000000e+00 4.110521e-01])

"""
function tang(psi::RV3)
    psi = psi.v
    PSI2 = psi'*psi
    PSI=sqrt(PSI2)
    PSI > 1.0e-2 ? f_1 = (cos(PSI)-1.0)/PSI2 : f_1 = -0.5 + PSI2/24.0 - PSI2^2/720.0
    PSI > 1.0e-4 ? f_2 = (1.0-sin(PSI)/PSI)/PSI^2 : f_2 = 1.0/6.0 - PSI2/120.0 + PSI2^2/5040.0
    tps = [0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0]
    T = eye(3) + f_1*tps + f_2*tps*tps 
    
    return Mat3(T)
end

"""
    tang(psi::RP3)
        function computing the tangent operator T::Mat3 from the 
        Euler-Rodrigues parameters p::RP3 

"""
function tang(p::RP3)


    P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
    tps = [0.0 -p[3] p[2]; p[3] 0.0 -p[1]; -p[2] p[1] 0.0]

    return T = Mat3(1.0/P_0*eye(3) +  (0.25/P_0*tps- 0.5*eye(3))*tps)

end
