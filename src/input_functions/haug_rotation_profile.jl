
"""
    haug_rotation_profile

    constant velocity  after raising period T.
"""
function haug_rotation_profile(itime,h,params)

    t = (itime-1)*h

    T = params[2]

    omega = params[1]

    tau=t/T

    tau <= 1.0 ? val = omega*T*(tau^2/2 +(cos(2*pi*tau)-1)/(2*pi)^2) : val = omega*T*(tau-1/2)


    return val

end
