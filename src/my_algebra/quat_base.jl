Base.getindex(r::Quat,i::Int) = r.q[i]
Base.:-(r::Quat) = Quat([r[1];-r.q[2:4]])

function Quat(a::RV3)
    PSI=sqrt(a[1]^2+a[2]^2+a[3]^2)
    if PSI > PREC
        COEF=sin(PSI/2.0)/PSI
        q=Quat([cos(PSI/2.0); a[1]*COEF; a[2]*COEF; a[3]*COEF])
    else
        q=Quat()
    end
    return q
end

function Quat()
    return Quat([1.0;  0.0;  0.0; 0.0])
end

function Quat(u::Quat,v::Quat)
    p_0=u.q[1]
    p=Vec3(u.q[2:4])
    r_0=v.q[1]
    r=Vec3(v.q[2:4])
    s_0=p_0*r_0-dotp(p,r)
    s=p_0*r+r_0*p+crossp(p,r)
    return Quat([s_0;s[1];s[2];s[3]])
end

function RV3(a::Quat)
    SINP2=sqrt(a[2]^2+a[3]^2+a[4]^2)
    if SINP2 > PREC
        COSP2=a[1]
        COEF=2.0*atan(SINP2,COSP2)/SINP2
    #    if  COEF > pi
    #        COEF -= 2.0*pi
    #    end
        b= RV3(COEF*a.q[2:4])
    else
        b = RV3(2.0*a.q[2:4])
    end
    return b
end
