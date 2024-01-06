
"""
    RP3

```math
\\underline{p} = 2.0 * \\underline{n} \\sin(\\phi/2)  =  RV3(p_1, p_2, p_3)
```
Mutable structure hosting a (3x1) set of Rodrigues - Euler parameters ``\\underline{p}``.

"""
mutable struct RP3
    v::Vector{Float64}
end


RP3() = RP3(zeros(3))

RP3(a::Vec3) = RP3(a.v)

RP3(x::Float64,y::Float64,z::Float64) =  RP3(Vector{Float64}(x,y,z))
