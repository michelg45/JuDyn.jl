function pull_vectors(loc::Array{Int,1},y::Vector{Float64})

    """
    function pull_vectors(loc::Array{Int,1},y::Vector)

        version 1: extraction of one vector

        input
            loc: localization vector
            y:   source vector

        output
            v = y[loc]

    """

        ndim = size(loc,1)
        v = Vector{Float64}(zeros(ndim))

        if findmin(loc)[1] > 0

            v = y[loc]
        
        else
        
            for i = 1:ndim
                loc[i] != 0 && (v[i] = y[loc[i]])
            end
            
        end
   
      
        return v
    end

function pull_vectors(loc::Array{Int,1},y1::Vector{Float64},y2::Vector{Float64})

      """
      function pull_vectors(loc::Array{Int,1},y::Vector,,y2::Vector)

          version 2: extraction of two vectors

          input
              loc:    localization vector
              y1, y2: source vectors

          output
              v1 = y1[loc], v2 = y2[loc]

      """

      ndim = size(loc,1)
      v1 = Vector{Float64}(zeros(ndim))
      v2 = Vector{Float64}(zeros(ndim))
    
      if findmin(loc)[1] > 0

        v1 = y1[loc]
        v2 = y2[loc]
    
      else
    
        for i = 1:ndim
            loc[i] != 0 && (v1[i] = y1[loc[i]]; v2[i] = y2[loc[i]])
        end
        
      end
     
      return (v1, v2)
  end

function pull_vectors(loc::Array{Int,1},y1::Vector{Float64},y2::Vector{Float64},y3::Vector{Float64})

  """
  function pull_vectors(loc::Array{Int,1},y::Vector,,y2::Vector,,y3::Vector)

      version 3: extraction of 3 vectors

      input
          loc:        localization vector
          y1, y2, y3: source vectors

      output
          v1 = y1[loc], v2 = y2[loc], v3 = y3[loc]

  """

  ndim = size(loc,1)
  v1 = Vector{Float64}(zeros(ndim))
  v2 = Vector{Float64}(zeros(ndim))
  v3 = Vector{Float64}(zeros(ndim))

  if findmin(loc)[1] > 0

    v1 = y1[loc]
    v2 = y2[loc]
    v3 = y3[loc]

  else

    for i = 1:ndim
        loc[i] != 0 && (v1[i] = y1[loc[i]]; v2[i] = y2[loc[i]] ; v3[i] = y3[loc[i]])
    end

  end



  return (v1, v2, v3)
end
