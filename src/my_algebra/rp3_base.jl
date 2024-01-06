Base.getindex(r::RP3,i::Int) = r.v[i]
Base.getindex(r::RP3,I) = [r.v[i] for i in I]
Base.:-(r::RP3) = RP3(-r.v)
Base.copy(r::RP3)= RP3(copy(r.v))

"""
    RP3(psi::RV3)
        function computing Euler-Rodrigues parameters p::RP3 from
        a Cartesian rotation vector psi::RV3

```math
        \\underline{p} = inv_p(\\underline{\\psi})
```       
"""
function RP3(psi::RV3)
    PSI2=psi[1]*psi[1]+psi[2]*psi[2]+psi[3]*psi[3]
    PSI=sqrt(PSI2)
    if PSI <= PREC
        return b = RP3(psi.v)
    else
        return b = RP3(2*sin(PSI/2.0)/PSI*psi.v)
    end
end

"""
    RV3(psi::RP3)
        function computing a Cartesian rotation vector psi::RV3 
            from Euler-Rodrigues parameters  p::RP3 
```math
\\underline{\\psi} = inv_\\psi(\\underline{p})
```       
"""
function RV3(p::RP3)
    P = sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3])
    if P <= PREC
        return psi = RV3(p.v)
    else
        return psi = RV3(2*asin(P/2.0)/P*p.v)
    end
end

"""
    RP3(R::Mat3)
        function computing the Euler-Rodrigues parameters from a rotation operator R::Mat3
```math
\\underline{p} = inv_p(\\underline{\\underline{R}})
```       

"""
function RP3(R::Mat3)

    coef = 1.0/sqrt(R[1,1]+R[2,2]+R[3,3]+1.0)

    return p = RP3( coef*(R[3,2] - R[2,3]), coef*(R[1,3] - R[3,1]), coef*(R[2,1] - R[1,2]) )

end

"""
    Quat(p::RP3)
        function computing a quaternion q::Quat from Euler-Rodrigues parameters  p::RP3
```math
\\underline{q} = inv_q(\\underline{p})
```       
"""
Quat(p::RP3) = Quat([sqrt(1- (a[1]*a[1]+a[2]*a[2]+a[3]*a[3])); Vector(0.5*p.v)])


RP3(r::Quat) = RP3(2.0*r.q[2:4])

"""
        function Quat(r::Quat,p::RP3)

        Performs the rotation composition R(q) = R1(r)*R2(p)
        suitable if norm(p) is small (rotation increment)
        input: Quaternion r::Quat and Euler - Rodrigues parameters  p::RP3
        output: quaternion q::Quat

"""

Quat(r::Quat,p::RP3) = Quat(r,Quat(p))

"""
        function RP3(r::Quat,q::Quat)

        Performs the rotation composition R(p) = R1(r)*R2(q)
        suitable if norm(p) is small (comutation of small relative rotation)
        input: Quaternions r::Quat and q::Quat
        output: Euler - Rodrigues parameters  p::RP3

"""

RP3(r::Quat,q::Quat) = RP3(Quat(r,q))

"""
        function RP3(psi1::RV3,psi2::RV3)

        Performs the rotation composition R(p) = R(psi1)*R(psi2)
        suitable if norm(p) is small (comutation of small relative rotation)
        input: Quaternions psi1::RV3 and psi2::RP3
        output: Euler - Rodrigues parameters  p::RP3

"""

RP3(psi1::RV3,psi2::RV3) =  RP3(Quat(Quat(psi1),Quat(psi2)))

Vec3(p::RP3) =  Vec3(p.v)
