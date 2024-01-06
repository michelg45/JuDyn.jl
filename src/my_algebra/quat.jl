"""
    Quat 

         mutable 4-component structure describing a quaternion q::Quat
 
    Example

```@example
        q = Quat()

        Quat([1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

```    

"""
mutable struct Quat

    q::Vector{Float64}

    function Quat(q)
         abs(norm(q) -1.0) > PREC ?  error("Quat input: q has not correct norm") : new(q)
    end
end
