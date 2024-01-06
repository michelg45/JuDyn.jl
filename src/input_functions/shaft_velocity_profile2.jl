
"""
    shaft_velocity_profile2
"""
function shaft_velocity_profile2(itime::Int,h::Float64,params::Vector{Float64})

    omega = params[1]
    A1 = params[2]
    A2 = params[3]
    t1 = params[4]
    t2 = params[5]
    t3 = params[6]
    t_in = params[7]
    time = itime*h
    (time  <= t_in) && (val = 0.0)
    (time > t_in && time <= t1) &&  (val = A1*omega*(1.0 - cos(pi*(time-t_in)/(t1-t_in)))/2.0)
    (time > t3) && (val = A2*omega)
    (time > t1 && time <= t2) && (val = A1*omega)
    (time > t2 && time <= t3) && (val =omega*(A1 + (A2-A1)/2.0*(1.0 - cos(pi*(time-t2)/(t3-t2)))))
    return val*time 

end
