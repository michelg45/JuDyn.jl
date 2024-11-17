"""
    RV3

        Mutable structure hosting a (3x1) Cartesian Rotation Vector.

        Examples
```@example

        phi1 = RV3(0.0, 2.0, -1.0)
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])

        v = [.5,1.0,2.0]
        3-element Vector{Float64}:
            5.000000e-01
            1.000000e+00
            2.000000e+00
        phi2 = RV3(v)
        RV3([5.000000e-01, 1.000000e+00, 2.000000e+00])

        phi3 = RV3(phi1,phi2)   # composition rule
        RV3([2.504197e+00, 1.202713e+00, 6.239186e-02])

```

"""
mutable struct RV3
    v::Vector{Float64}
end

RV3() = RV3(zeros(3))

RV3(a::Vec3) = RV3(a.v)

RV3(x::Float64,y::Float64,z::Float64) =  RV3([x; y; z])
