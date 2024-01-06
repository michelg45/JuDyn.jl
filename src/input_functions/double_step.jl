"""
    double_step
"""
function double_step(itime,h,params)

    t = (itime-1)*h

    amp = params[1]
    t1 = params[2]
    t2 = params[3]
    t3 = params[4]

    val=0
        if (t > 0 && t <= t1)
            val=1
        elseif (t > t2 && t <= t3)
            val=-1
        end
    return val*amp
end
