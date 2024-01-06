"""
    linear_velocity
"""
function linear_velocity(itime::Int,h::Float64,params::Vector{Float64})

    acc = params[1]
    v_max = params[2]
    t1 = params[3]
    time = itime*h
    time < t1 ? omega = 0 : omega = acc*(time-t1)
    
    omega < v_max ? val = omega*(time-t1)/2.0 : (T = v_max/acc; val = 0.5*acc*T^2 + v_max*(time -T - t1))
    return val 

end
