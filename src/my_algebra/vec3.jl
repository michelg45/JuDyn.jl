"""
    Vec3

        mutable structure containing a 3-D vector (position, velocity, displacement)

        Examples 

```@example
        a = Vec3()
        Vec3([0.000000e+00, 0.000000e+00, 0.000000e+00])

        a = Vec3(2.0,4.0,5.0)
        Vec3([2.000000e+00, 4.000000e+00, 5.000000e+00])

        b = a.v
        3-element Vector{Float64}:
            2.000000e+00
            4.000000e+00
            5.000000e+00

        c = Vec3(b)
        Vec3([2.000000e+00, 4.000000e+00, 5.000000e+00]) 
```  
"""
mutable struct Vec3
    v::Vector{Float64}
end

function Vec3()
    return Vec3(zeros(3))
end

Vec3(x::Float64,y::Float64,z::Float64) =  Vec3([x; y; z])

