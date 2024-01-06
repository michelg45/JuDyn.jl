"""
    CurrentFrame

        Data structure hosting the current node reference frames during the incrementation process. 
"""
mutable struct CurrentFrame
    nbr::Int
    type::String
    inv_loc::Vector{Int}
    x::Vec3
    xn::Vec3
    Dx::Vec3
    xdot::Vec3
    psi::RV3
    psin::RV3
    psidot::Vec3
    Dpsi::RV3
    W::Vec3
    R::Mat3
    T::Mat3
    DW::Mat3

    function CurrentFrame()

        nbr=0
        type="frame"
        inv_loc=Vector{Int}(zeros(6))
        x=Vec3()
        xn=Vec3()
        Dx=Vec3()
        xdot=Vec3()
        psi=RV3()
        psin=RV3()
        Dpsi=RV3()
        psidot=Vec3()
        W=Vec3()
        R=Mat3()
        T=Mat3()
        DW=Mat3()

        new(nbr,type,inv_loc,x,xn,Dx,xdot,psi,psin,psidot,Dpsi,W,R,T,DW)

    end
end # CurrentFrame