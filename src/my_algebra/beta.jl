"""
    beta(x::Float64)
        function computing the series expansion of the function
\$\\beta(x) = 2.0*(1-\\cos(x))/x^2\$
"""
function beta(x::Float64)
  add = 0.5
  i = 0
  f = add
  while abs(add) >= PRC
      i += 1
      add = -add*x^2/((2*i+1)*(2*i+2))
      f += add
   end
  return 2.0*f
end
