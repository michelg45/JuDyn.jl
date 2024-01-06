"""
    shaft_velocity_profile
"""
function shaft_velocity_profile(itime::Int,h::Float64,params::Vector{Float64})

    omega = params[1]
    A1 = params[2]
    A2 = params[3]
    t1 = params[4]
    t2 = params[5]
    t3 = params[6]
    time = itime*h
    v1 = -(A1*t1*omega)/2.0
    v2 = A1*t2*omega + v1
    v3 = -(((A2-A1)*t2)/2.0+A1*t2)*omega +v2
    v4 = (((A2-A1)*t3)/2.0+A1*t3)*omega + v3 - A2*t3*omega
    (time <= t1) &&  (val = (A1*(time-(sin((pi*time)/t1)*t1)/pi)*omega)/2.0)
    (time > t1 && time <= t2) && (val =  A1*time*omega +v1) 
    (time > t2 && time <= t3) && (val =  (((A2-A1)*(time-((t3-t2)*sin((pi*(time-t2))/(t3-t2)))/pi))/2.0 + A1*time)*omega + v3) 
    (time > t3) && (val = A2*time*omega + v4)
        
    return val 

end
