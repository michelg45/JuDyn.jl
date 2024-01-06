
"""
    element_sparse_indexes
"""
function element_sparse_indexes(inv_loc::Int)
    index_I = Vector{Int}[]
    index_J = Vector{Int}[]
    for i = 1:N
       index_I =  [index_I; [inv_loc[i] for  j in 1:N]]
       index_J =  [index_J; inv_loc]
    end
    return index_I, index_J
end
    