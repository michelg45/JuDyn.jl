"""
    DinvtangT(a::Vec3,psi::RV3)

        function computing the D::Mat3 matrix

```math
\\underline{\\underline{D}}(\\underline{a},\\underline{\\psi}) = D_\\psi (T(\\underline{\\psi})^{-T}\\underline{a})
``` 

"""
function DinvtangT(a::Vec3,psi::RV3)

    PSI2=psi[1]*psi[1]+psi[2]*psi[2]+psi[3]*psi[3]
    PSI=sqrt(PSI2)
    T = Mat3(0.5*[0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0])
    if  PSI > PREC
        bracket = a.v*transpose(psi.v)+psi.v*transpose(a.v)-2*transpose(psi.v)*a.v*eye(3)
        delta = (1 - PSI/2*cot(PSI/2))/PSI2
        T +=  Mat3(delta*bracket + (delta^2 + 3*delta/PSI2 + PSI2/4.0)*(crossp(psi,crossp(psi,a)).v*transpose(psi.v)))
    end
    return T
end

"""
    Dinvtang(a::Vec3,psi::RV3)

        function computing the D::Mat3 matrix

```math
\\underline{\\underline{D}}(\\underline{a},\\underline{\\psi}) = D_\\psi (\\underline{\\underline{T}}(\\underline{\\psi})^{-1}\\underline{a})
``` 

"""
function Dinvtang(a::Vec3,psi::RV3)

    PSI2=psi[1]*psi[1]+psi[2]*psi[2]+psi[3]*psi[3]
    PSI=sqrt(PSI2)
    T = Mat3(0.5*[0.0 a[3] -a[2]; -a[3] 0.0 a[1]; a[2] -a[1] 0.0])
    if  PSI > PREC
        bracket = a.v*transpose(psi.v)+psi.v*transpose(a.v)-2*transpose(psi.v)*a.v*eye(3)
        delta = (1 - PSI/2*cot(PSI/2))/PSI2
        T +=  Mat3(delta*bracket + (delta^2 + 3*delta/PSI2 + PSI2/4.0)*(crossp(psi,crossp(psi,a)).v*transpose(psi.v)))
    end
    return T
end

"""
    Dinvtang(a::Vec3,p::RP3)

        function computing the D::Mat3 matrix

```math
\\underline{\\underline{D}}(\\underline{a},\\underline{p}) = D_\\psi (\\underline{\\underline{T}}(\\underline{p})^{-1}\\underline{a})
``` 

"""
function Dinvtang(u::Vec3, p::RP3)
    P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
       return 0.5*tilde(u) -1/(4.0*P_0)*Mat3(u,Vec3(p))
end
