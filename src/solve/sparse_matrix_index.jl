"""
    sparse_index
"""
function sparse_index(loc1::Vector{Int}, loc2::Vector{Int})

    n1 = size(loc1,1)
    n2 = size(loc2,1)

    I = Vector{Int64}[]
    J = Vector{Int64}[]

    for j = 1:n2
        I  = Vector{Int64}([I; loc1])
        J  = Vector{Int64}([J; [loc2[j] for i in 1:n1]])
    end

return I, J

end
