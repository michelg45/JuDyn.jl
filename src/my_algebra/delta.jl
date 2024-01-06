"""
    delta(x::Float64)   
        function computing the series expansion of
\$\\delta(x) = (1-\\gamma(x))/x^2\$
"""
function delta(x::Float64)
 i = 1
 f = 1.0/12.0
 add = 0.5
 xpower = 1.
 while abs(add) >= PRC && i < 15
     i += 1
     xpower = -xpower*x^2
     add =xpower*coef_gamma(i)
     f += add
 end
 return f
end
