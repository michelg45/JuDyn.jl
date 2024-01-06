
"""
    inverse_step
"""
function inverse_step(itime,h,params)

    t = (itime-1)*h

    amp = params[1]
    t1 = params[2]

    t <= t1 ? val= params[1] : val=0.

    return val
end