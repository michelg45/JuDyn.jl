function vector_mean(v)
      sum = 0.
      N = size(v,1)
      if size(v,2) != 1
            error(println("input is not a vector"))
      end
      for i=1:N
            sum +=v[i]
      end
      return sum/N
end
