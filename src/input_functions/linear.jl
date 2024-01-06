
"""
    linear
"""
function linear(itime::Int,h::Float64,params::Vector{Float64})

    T = params[1]
    amp = params[2]
    return val = amp*(itime-1)*h/T

end
