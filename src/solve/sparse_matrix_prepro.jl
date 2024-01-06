"""
    sparse_matrix_prepro
"""
function sparse_matrix_prepro()

    mc = Main.model_container
    ec = Main.element_container

    Nel = mc.Elements
    Ndofs_q = mc.Ndofs_x + mc.Ndofs_int + mc.Ndofs_mult
    Ndofs_v = mc.Ndofs_v

    rotation = mc.uniform_rotation

    iptr_elem_qq = zeros(Int,Nel)
    iptr_elem_vq = zeros(Int,Nel)
    iptr_elem_qv = zeros(Int,Nel)

    I_index_qq = Vector{Int}[]
    J_index_qq = Vector{Int}[]
    subset_q = Vector{Int}[]
    I_index_vq = Vector{Int}[]
    J_index_vq = Vector{Int}[]
    I_index_qv = Vector{Int}[]
    J_index_qv = Vector{Int}[]


    ip_qq = 1
    ip_vq = 1
    ip_qv = 1


    for iel =1:Nel

        if ec.element_types[iel] != "node_force"
            inv_loc_x = ec.inv_loc_x[iel]
            inv_loc_int = ec.inv_loc_int[iel]
            inv_loc_mult = ec.inv_loc_mult[iel]
            inv_loc_v = nonnuls(ec.inv_loc_v[iel])
            inv_loc_q = nonnuls([inv_loc_x;inv_loc_int;inv_loc_mult])

            # println("element ", iel, "inv_loc_q ", inv_loc_q)
            # println("element ", iel, "inv_loc_v ", inv_loc_v)

            I, J = sparse_index(inv_loc_q,inv_loc_q)
            I_index_qq = [I_index_qq; I]
            J_index_qq = [J_index_qq; J]
            iptr_elem_qq[iel] = ip_qq
            length = size(I,1)
            ip_qq += length

            if size(inv_loc_v,1) > 0

                I, J = sparse_index(inv_loc_v,inv_loc_q)
                I_index_vq = [I_index_vq; I]
                J_index_vq = [J_index_vq; J]
                iptr_elem_vq[iel] = ip_vq
                length = size(I,1)
                subset_q = [subset_q; inv_loc_q]
                ip_vq += length

            end

        end

    end

    nz_qq = size(J_index_qq,1)
    nz_vq = size(J_index_vq,1)

    if nz_vq > 0
        subset_q = unique(x -> x, subset_q)
    end

    global S = SparseMatrixStructure(Nel,Ndofs_q,Ndofs_v,nz_qq,nz_vq)

    S.I_index_qq = I_index_qq
    S.J_index_qq = J_index_qq
    S.iptr_elem_qq = iptr_elem_qq
    S.subset_q = subset_q

    if Ndofs_v > 0
        S.I_index_vq = I_index_vq
        S.J_index_vq = J_index_vq
        S.iptr_elem_vq = iptr_elem_vq
    end


    println("sparse matrix structure : ")
    println("number of kinematic dofs ", Ndofs_q)
    println("number of velocity dofs ", Ndofs_v)
    println("Vector S_qq size ", nz_qq)
    println("Vector S_vq size ", nz_vq)
    println("subset_q size ", size(subset_q,1))
    # println("Vector S_qq pointers ", iptr_elem_qq)
    # println("Vector S_vq pointers ", iptr_elem_vq)


    return S

end
