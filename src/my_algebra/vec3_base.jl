Base.copy(r::Vec3)=Vec3(copy(r.v))
Base.getindex(r::Vec3,i::Int) = r.v[i]
Base.getindex(r::Vec3,I) = [r.v[i] for i in I]
Base.setindex!(a::Vec3,val,i::Int) =  a.v[i] = val
Base.:-(a::Vec3) = Vec3(-a.v)
Base.:+(a::Vec3, b::Vec3) = Vec3(a.v+b.v)
Base.:-(a::Vec3, b::Vec3) = Vec3(a.v-b.v)
Base.:*(a::Float64, b::Vec3) =Vec3(a*b.v)
Base.:*(a::Int, b::Vec3) =Vec3(a*b.v)
Base.:*(b::Vec3,a::Float64) =Vec3(a*b.v)
Base.:/(a::Vec3,b::Float64) =Vec3(a.v/b)
Base.:*(a::Vec3, b::Vec3) = a.v[1]*b.v[1] + a.v[2]*b.v[2]+a.v[3]*b.v[3]

"""
    norm2(a::Vec3)

    function calculating the norn of a vector of type ::Vec3.

       calling sequence: nr = norm2(a) 

       Examples
    
       a = Vec3(5.0, 2.0, -1.0)
        Vec3([5.000000e+00, 2.000000e+00, -1.000000e+00])
       nr = norm2(a)
       5.477226e+00

       b = a.v
       3-element Vector{Float64}:
        5.000000e+00
        2.000000e+00
        -1.000000e+00
       nr = norm2(b)
         ERROR: MethodError: no method matching norm2(::Vector{Float64})
         Closest candidates are: 
         norm2(::Vec3)
    
"""
norm2(a::Vec3)= sqrt(a.v[1]*a.v[1] + a.v[2]*a.v[2]+a.v[3]*a.v[3])

"""
    dotp(a::Vec3,b::Vec3)

    function calculating the dot product of two vectors of type ::Vec3.

       calling sequence: c = dotp(a,b) 

       Example
    
        a = Vec3(5.0, 2.0, -1.0)
        Vec3([5.000000e+00, 2.000000e+00, -1.000000e+00])
        b = Vec3(2.0,4.0,-2.0)
        Vec3([2.000000e+00, 4.000000e+00, -2.000000e+00])
        c = dotp(a,b)
        2.000000e+01
  
"""
dotp(a::Vec3, b::Vec3) = a.v[1]*b.v[1] + a.v[2]*b.v[2]+a.v[3]*b.v[3]

"""
    crossp(a::Vec3,b::Vec3)

    function calculating the crossp product of two vectors of type ::Vec3.

       calling sequence: c = crossp(a,b) 

       Example
    
        a = Vec3(5.0, 2.0, -1.0)
        Vec3([5.000000e+00, 2.000000e+00, -1.000000e+00])
        b = Vec3(2.0,4.0,-2.0)
        Vec3([2.000000e+00, 4.000000e+00, -2.000000e+00])
        c = crossp(a,b)
        Vec3([0.000000e+00, 8.000000e+00, 1.600000e+01])
  
"""
crossp(a::Vec3, b::Vec3) =  Vec3(-a.v[3]*b.v[2]+a.v[2]*b.v[3],-a.v[1]*b.v[3]+a.v[3]*b.v[1], -a.v[2]*b.v[1]+a.v[1]*b.v[2]) 
