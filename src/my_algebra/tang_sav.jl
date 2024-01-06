"""
    tang(psi::RV3,a::Vec3)
        function computing the product of a vector a::Vec3 
        by a tangent operator with Cartesian vector psi::RV3
```math
    \\underline{b} = \\underline{\\underline{T}}(\\underline{\\psi})  \\underline{a} 
```
"""
function tang(psi::RV3,a::Vec3)
    psi = psi.v
    a = a.v
    PSI=sqrt(psi'*psi)
    PSI <= PREC ?  b = a - 0.5*cross(psi,a) : b = sin(PSI)/PSI*a + (cos(PSI)-1)/PSI^2*cross(psi,a) + (1.0-sin(PSI)/PSI)/PSI^2*dot(psi,a)*psi
    return Vec3(b)
end

"""
    tang(psi::RV3)
        function computing the tangent operator T::Mat3
        from the Cartesian rotation vector psi::RV3 
```math
\\underline{\\underline{T}} = \\underline{\\underline{T}}(\\underline{\\psi}) 
```
"""
function tang(psi::RV3)
    psi = psi.v
    PSI= norm(psi)
    tps = [0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0]
    PSI <= PREC ? T = eye(3)-0.5*tps : T = sin(PSI)/PSI*eye(3) + (cos(PSI)-1)/PSI^2*tps + ((1.0-sin(PSI)/PSI)/PSI^2*psi)*psi'
    
    return Mat3(T)
end

"""
    tang(psi::RP3)
        function computing the tangent operator T::Mat3 from the 
        Euler-Rodrigues parameters p::RP3 
```math
\\underline{\\underline{T}} = \\underline{\\underline{T}}(\\underline{p}) 
```
"""
function tang(p::RP3)


    P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
    tps = [0.0 -p[3] p[2]; p[3] 0.0 -p[1]; -p[2] p[1] 0.0]

    return T = Mat3(1.0/P_0*eye(3) +  (0.25/P_0*tps- 0.5*eye(3))*tps)

end
