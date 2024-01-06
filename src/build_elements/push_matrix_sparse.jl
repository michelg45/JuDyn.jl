"""
    push_matrix_sparse
"""
function push_matrix_sparse(iel::Int,loc_q::Array{Int,1},loc_v::Array{Int,1},S_el::Array{Float64,2})

    nz_v = size(nonnuls(loc_v),1)
    if nz_v > 0
        sps = Main.Solve.sparse_matrix_S
        N_q = size(loc_q,1)
        N_v = size(loc_v,1)
        locel_q = nonnuls((loc_q .> 0).*[i for i in 1:N_q])
        locel_v = nonnuls((loc_v .> 0).*[i+N_q for i in 1:N_v])
        loc_q = nonnuls(loc_q)
        loc_v = nonnuls(loc_v) .+ sps.Ndofs_q


        length_qq  = size(locel_q,1)^2
        ip_qq = sps.iptr_elem_qq[iel]
        range_qq = [ip_qq+i-1 for i in 1:length_qq]

        sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])

        length_vq  = size(locel_q,1)*size(locel_v,1)
        ip_vq = sps.iptr_elem_vq[iel]
        range_vq = [ip_vq+i-1 for i in 1:length_vq]
        sps.S_vq[range_vq] = vec(S_el[locel_v,locel_q])
        sps.S_qvt[range_vq] = vec(transpose(S_el[locel_q,locel_v]))
        
    else
        push_matrix_sparse(iel,loc_q,S_el)
    end

    return

end

function push_matrix_sparse(iel::Int,loc_q::Array{Int,1},S_el::Array{Float64,2})

    sps = Main.Solve.sparse_matrix_S
    N_q = size(loc_q,1)
    locel_q = nonnuls((loc_q .> 0).*[i for i in 1:N_q])
    loc_q = nonnuls(loc_q)
    length_qq  = size(locel_q,1)^2
    ip_qq = sps.iptr_elem_qq[iel]
    range_qq = [ip_qq+i-1 for i in 1:length_qq]
    sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])

    return

end