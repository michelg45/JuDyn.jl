"""
    alpha(x::Float64)
        function computing the series expansion of the function
\$\\alpha(x) = \\sin(x)/x\$
"""
function alpha(x::Float64)
    add = 1.0
    i = 0
    f = 1.0
    while abs(add) >= PRC
        i += 1
        add = -add*x^2/(2*i*(2*i+1))
        f += add
    end
    return f
end
