"""
    condensed_element_matrix    
"""
function condensed_element_matrix(iel::Int,inv_loc_q::Array{Int,1},inv_loc_v::Array{Int,1},S_el::Array{Float64,2})
    
    nonnuls(a::Vector{Int}) = a[a .> 0]
    nz_v = size(nonnuls(inv_loc_v),1)
    if nz_v > 0
        N_q = size(inv_loc_q,1)
        N_v = size(inv_loc_v,1)
        locel_q = nonnuls((inv_loc_q .> 0).*[i for i in 1:N_q])
        locel_v = nonnuls((inv_loc_v .> 0).*[i+N_q for i in 1:N_v])
        inv_loc =  [nonnuls(inv_loc_q); nonnuls(inv_loc_v) .+ Main.model_container.Ndofs_q]
        locel = [locel_q; locel_v]
    else
        N_q = size(inv_loc_q,1)
        locel = nonnuls((inv_loc_q .> 0).*[i for i in 1:N_q])
        inv_loc =  nonnuls(inv_loc_q)
    end

    S_el = S_el[locel,locel]
        
    return inv_loc, S_el

end

function condensed_element_matrix(iel::Int,inv_loc::Array{Int,1},S_el::Array{Float64,2})
    
    N = size(inv_loc,1)
    locel = nonnuls((inv_loc .> 0).*[i for i in 1:N])
    inv_loc =  nonnuls(inv_loc)

    S_el = S_el[locel,locel]
        
    return inv_loc, S_el

end