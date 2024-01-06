"""
    push_element_sparse

            function to push element entities (element array 'S_el', residual vector 'res_el', load vector 'p_el')
            into the corresponding structural sets 'S', 'res' and 'p'. The structural iteration matrix S is stored 
            in sparse form in the sparse array "Main.Solve.sparse_matrix_S".

            calling sequences:

                push_element_sparse(res,p,iel,loc_q,loc_v,S_el, res_el,p_el)

                push_element_sparse(res,p,iel,loc_q,loc_v,S_el, res_el)

                push_element_sparse(res,p,iel,loc_q,S_el, res_el,p_el)

                push_element_sparse(res,p,iel,loc_q,S_el, res_el)


"""
function push_element_sparse(res::Vector{Float64},p::Vector{Float64},iel::Int,loc_q::Array{Int,1},loc_v::Array{Int,1},S_el::Array{Float64,2},res_el::Vector{Float64}, p_el::Vector{Float64})

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

        res[loc_q]   += (res_el[locel_q]+ S_el[locel_q,locel_v]*res_el[locel_v])

        S_el[locel_q,locel_q] += S_el[locel_q,locel_v]*S_el[locel_v,locel_q]
        sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])

        res[loc_v]   =  res_el[locel_v]
        p[loc_q]     += p_el[locel_q]
        length_vq  = size(locel_q,1)*size(locel_v,1)
        ip_vq = sps.iptr_elem_vq[iel]
        range_vq = [ip_vq+i-1 for i in 1:length_vq]
        sps.S_vq[range_vq] = vec(S_el[locel_v,locel_q])
    else
        push_element_sparse(res,p,iel,loc_q,S_el,res_el, p_el)
    end

    return

end

function push_element_sparse(res::Vector{Float64},iel::Int,loc_q::Array{Int,1},loc_v::Array{Int,1},S_el::Array{Float64,2},res_el::Vector{Float64})

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

        res[loc_q]   += (res_el[locel_q]+ S_el[locel_q,locel_v]*res_el[locel_v])

        S_el[locel_q,locel_q] += S_el[locel_q,locel_v]*S_el[locel_v,locel_q]
        sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])

        res[loc_v]   =  res_el[locel_v]
        length_vq  = size(locel_q,1)*size(locel_v,1)
        ip_vq = sps.iptr_elem_vq[iel]
        range_vq = [ip_vq+i-1 for i in 1:length_vq]
        sps.S_vq[range_vq] = vec(S_el[locel_v,locel_q])
    else
        push_element_sparse(res,iel,loc_q,S_el,res_el)
    end
    return

end

function push_element_sparse(res::Vector{Float64},p::Vector{Float64},iel::Int,loc_q::Array{Int,1},S_el::Array{Float64,2},res_el::Vector{Float64}, p_el::Vector{Float64})

    sps = Main.Solve.sparse_matrix_S
    N_q = size(loc_q,1)
    locel_q = nonnuls((loc_q .> 0).*[i for i in 1:N_q])
    loc_q = nonnuls(loc_q)
    length_qq  = size(locel_q,1)^2
    ip_qq = sps.iptr_elem_qq[iel]
    range_qq = [ip_qq+i-1 for i in 1:length_qq]
    sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])
    res[loc_q]   += res_el[locel_q]
    p[loc_q]     += p_el[locel_q]


    return

end

function push_element_sparse(res::Vector{Float64},iel::Int,loc_q::Array{Int,1},S_el::Array{Float64,2},res_el::Vector{Float64})

    sps = Main.Solve.sparse_matrix_S
    N_q = size(loc_q,1)
    locel_q = nonnuls((loc_q .> 0).*[i for i in 1:N_q])
    loc_q = nonnuls(loc_q)
    length_qq  = size(locel_q,1)^2
    ip_qq = sps.iptr_elem_qq[iel]
    range_qq = [ip_qq+i-1 for i in 1:length_qq]
    sps.S_qq[range_qq] = vec(S_el[locel_q,locel_q])
    res[loc_q]   += res_el[locel_q]

    return

end
