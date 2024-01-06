"""
    coef_gamma(n::Int)
        function computing the nth coefficient of the series expansion of
                    gamma(x) = x/2*cot(x/2)
        could be used to build a constant table for computing gamma(x)
        and associated functions.
"""
function coef_gamma(n::Int)

      factor = 1.
      for i=2:2*n
          factor = factor*i
      end
      B = bernoulli(2*n)
      num = Float64(numerator(B))
      den = Float64(denominator(B))*factor
      return num/den
  end

"""
    gamma(x::Float64)
        function computing the series expansion of
\$\\gamma(x) = x/2*\\cot(x/2)\$
"""
function gamma(x::Float64)
   i = 0
   f = 1.
   add = 1.
   xpower = 1.0
   fact = 1.0
   while abs(add) >= PRC && i <= 15
       i += 1
       xpower = -xpower*x^2
       add = xpower*coef_gamma(i)
       f += add
    end
   return f
end
