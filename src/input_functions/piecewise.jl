"""
    piecewise
"""
function piecewise(itime::Int,h::Float64,params)

    N_params = size(params,1)
    mod(N_params,2) != 0 &&  error("number of parameters of piecewise function must be even")
    N_steps = Int(N_params/2)
    i_steps = [(i-1)*2+1 for i in 1:N_steps] 
    i_vals = i_steps .+ 1
    t = params[i_steps]
    vals = params[i_vals]

    tim = (itime-1)*h

    val = 0.

    for i = 2:N_steps
        if  tim >= t[i-1] &&  tim <= t[i]
            a = (vals[i] - vals[i-1])/(t[i] - t[i-1])
            b = vals[i] -a*t[i]
            val = a*tim + b
        end
    end

    return val

end
function piecewise(tim::Float64,params)

    N_params = size(params,1)
    mod(N_params,2) != 0 &&  error("number of parameters of piecewise function must be even")
    N_steps = Int(N_params/2)
    i_steps = [(i-1)*2+1 for i in 1:N_steps] 
    i_vals = i_steps .+ 1
    t = params[i_steps]
    vals = params[i_vals]

       val = 0.

    for i = 2:N_steps
        if  tim >= t[i-1] &&  tim <= t[i]
            a = (vals[i] - vals[i-1])/(t[i] - t[i-1])
            b = vals[i] -a*t[i]
            val = a*tim + b
        end
    end

    return val

end