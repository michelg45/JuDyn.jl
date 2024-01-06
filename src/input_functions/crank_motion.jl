
"""
    crank_motion
"""
function crank_motion(itime,h,params)

    # torque  applied a node 1 (hub)

    t = (itime-1)*h

    T = params[1] # period

    tau=t/T
    if ( tau <= 1)
        val= pi*(1-cos(pi*tau))/2.0
    end
    if ( tau > 1)
        val= 1.0*pi
    end


return val

end
