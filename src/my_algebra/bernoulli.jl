"""
    bernoulli2(n::Int)
        function computing the Float64 value of the Bernoulli number B(n)
        needed for the series expansion of gamma(x).
"""
function bernoulli2(n::Int)

      A = Vector{Float64}(undef, n + 1)
      for m = 0 : n
          A[m + 1] = 1 / (m + 1)
          for j = m : -1 : 1
              A[j] = j * (A[j] - A[j + 1])
          end
      end
      return A[1]
  end

"""
    bernoulli(n::Int)
        function computing the rational expression of the Bernoulli number B(n)
        needed for the series expansion of gamma(x).
        num = numerator(A[1])
        den = denominator(A[1])
"""
  function bernoulli(n::Int)

      A = Vector{Rational{BigInt}}(undef, n + 1)
      for m = 0 : n
          A[m + 1] = 1 // (m + 1)
          for j = m : -1 : 1
              A[j] = j * (A[j] - A[j + 1])
          end
      end
      return A[1]
  end
