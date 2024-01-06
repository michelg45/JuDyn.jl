"""
    RV3

        Mutable structure hosting a (3x1) Cartesian Rotation Vector.

        Example

        phi = RV3(0.0, 2.0, -1.0)
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])


"""
mutable struct RV3
    v::Vector{Float64}
end

RV3() = RV3(zeros(3))

RV3(a::Vec3) = RV3(a.v)

RV3(x::Float64,y::Float64,z::Float64) =  RV3([x; y; z])
