"""
    Quat 

         mutable 4-component structure describing a quaternion 
 
    Examples

```@example
        q1 = Quat()
        Quat([1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

        phi = RV3(0.0, 2.0, -1.0)
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])
        q2 = Quat(phi)
        Quat([4.374512e-01, 0.000000e+00, 8.043066e-01, -4.021533e-01])

        q3 = Quat(q1,q2)    # composition rule
        Quat([4.374512e-01, 0.000000e+00, 8.043066e-01, -4.021533e-01])
```   

"""
mutable struct Quat

    q::Vector{Float64}

    function Quat(q)
         abs(norm(q) -1.0) > PREC ?  error("Quat input: q has not correct norm") : new(q)
    end
end
