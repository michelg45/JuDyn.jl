"""
    Dtang(a::Vec3,b::RV3)
        function computing the D::Mat3 matrix
```math
\\underline{\\underline{D}}(\\underline{a},\\underline{b}) = D_b(\\underline{\\underline{T}}(\\underline{b})\\underline{a})
``` 
        It is currently approximated as 
```math
        0.5 \\; \\tilde{\\underline{a}}
```
"""
function Dtang(a::Vec3,b::RV3)

    DbTba = 0.5*tilde(a)

   """ 
    b = Vec3(b.v)
    b2 = b[1]^2 + b[2]^2 + b[3]^2
    nb = sqrt(b2)

    nb > PREC && (
        tb = tilde(b);
        alpha = sin(nb)/nb; 
        beta = 2.0*(1-cos(nb))/b2;
        gamma = nb/2.0*cot(nb/2.0);
        ab = transpose(a.v)*b.v;
        DbTba  += (1.0 - gamma)/b2*(outerp(a,b) + outerp(b,a) - 2.0*ab*I3)  + (1/beta + gamma -2.0)/b2^2*ab*tb*tb
    )
     
    """
    return DbTba
end

    
