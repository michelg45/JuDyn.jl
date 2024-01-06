"""
    triangle_impulse
"""
function triangle_impulse(itime,h,params)
    ampl = params[1]
    t1 = params[2]
    t2 = params[3]
    t3 = params[4]

    time = itime*h

        val = 0.
    if  time > t1 &&  time <= t2
        val = (time-t1)/(t2-t1)
    elseif  time > t2 &&  time < t3
        val   =  (t3 - time)/(t3-t2)
    end
    return val*ampl
end
