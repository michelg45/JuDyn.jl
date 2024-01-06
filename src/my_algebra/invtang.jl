"""
    invtang(psi::RV3)
    
        function computing the inverse of the tangent operator T(psi).

        calling sequence: invT = invtang(psi)

        input: rotation vector psi::RV3
        output: inverse tangent operator T(psi)^{-1}.
        PREC: precision set to machine precision eps**1/2.
        if (PSI < PREC)  T is limited to its linearized expression.

"""
function invtang(psi::RV3)

    psi = psi.v
    PSI2 = psi'*psi
    PSI = sqrt(PSI2)
    tps=[0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0]
    PSI > 1.e-4 ? f_3 = (1.0 - PSI/2.0*cot(PSI/2.0))/PSI2 : f_3 = 1.0/12.0 + PSI2/720.0 + PSI2*PSI2/30240.0

    return Mat3(eye(3) + 0.5*tps + f_3*tps*tps)

end

"""
    invtang(psi::RV3,a::Vec3)

        function computing the product T^{-1}a of a vector a::Vec3 by the inverse of the tangent operator T(psi).

        calling sequence: b = invtang(psi,a)

        input: rotation vector psi::RV3 and vector a::Vec3.
        output: vector b::Vec3.

        PREC: precision set to machine precision eps**1/2.
        if (PSI < PREC)  T is limited to its linearized expression.

"""
function invtang(psi::RV3,a::Vec3)

    psi = psi.v
    a = a.v
    PSI2 = psi'*psi
    PSI = sqrt(PSI2)

    PSI > 1.e-2 ? f_3 = (1.0 - PSI/2.0*cot(PSI/2.0))/PSI2 : f_3 = 1.0/12.0 + PSI2/720.0 + PSI2*PSI2/30240.0

    return Vec3(a + 0.5*cross(psi,a) + f_3*cross(psi,cross(psi,a)))

end

"""
    invtang(p::RP3)

    function computing the inverse of the tangent operator T(p).

        calling sequence: invT = invtang(p)

        input: Euler - Rodrigues parameters  p::RP3
        output: inverse of tangent operator T(p)^{-1}.

"""
function invtang(p::RP3)


        P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
        tps = [0.0 -p[3] p[2]; p[3] 0.0 -p[1]; -p[2] p[1] 0.0]

        return T = Mat3(P_0*eye(3) + 0.5*tps)

end
