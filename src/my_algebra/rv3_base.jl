

Base.getindex(r::RV3,i::Int) = r.v[i]
Base.getindex(r::RV3,I) = [r.v[i] for i in I]
Base.setindex!(r::RV3,b::Float64,i::Int) = (r.v[i]=b)
Base.:-(a::RV3) = RV3(-a.v)
Base.copy(r::RV3)= RV3(copy(r.v))

"""
    crossp(a::RV3,b::Vec3)

    function calculating the crossp product of a vector of type ::RV3 with a vector of type ::Vec3.

       calling sequence: c = crossp(a,b) 

       Example
    
        a = RV3(5.0, 2.0, -1.0)
        RV3([5.000000e+00, 2.000000e+00, -1.000000e+00])
        b = Vec3(2.0,4.0,-2.0)
        Vec3([2.000000e+00, 4.000000e+00, -2.000000e+00])
        c = crossp(a,b)
        Vec3([0.000000e+00, 8.000000e+00, 1.600000e+01])

        d = crossp(b,a)
        ERROR: MethodError: no method matching crossp(::Vec3, ::RV3)
        Closest candidates are:
        crossp(::Vec3, ::Vec3)

"""
crossp(u::RV3,v::Vec3) = Vec3(v[3]*u[2]-v[2]*u[3], v[1]*u[3]-v[3]*u[1] ,v[2]*u[1]-v[1]*u[2])

"""
    tilde(psi::RV3)

    function constructing the skew-symmetric matrix (type ::Mat3) associated to a rotation vector (type ::RV3).

       calling sequence: A = tilde(psi)
       
       example:

       psi = RV3(5.0, 2.0, -1.0)
       A = tilde(psi)
       Mat3([0.000000e+00 1.000000e+00 2.000000e+00; -1.000000e+00 0.000000e+00 -5.000000e+00; 
       -2.000000e+00 5.000000e+00 0.000000e+00])

"""
tilde(psi::RV3) = Mat3([0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0])
Base.:*(a::Float64,r::RV3) = RV3(a*(r.v))