"""
    SparseMatrixStructure
"""
mutable struct SparseMatrixStructure
    Ndofs_q::Int
    Ndofs_v::Int
    I_index_qq::Vector{Int}
    J_index_qq::Vector{Int}
    I_index_vq::Vector{Int}
    J_index_vq::Vector{Int}
    iptr_elem_qq::Vector{Int}
    iptr_elem_vq::Vector{Int}
    subset_q::Vector{Int}
    S_qq::Vector{Float64}
    S_vq::Vector{Float64}
    S_qvt::Vector{Float64}
 
    function SparseMatrixStructure(Nel,Ndofs_q,Ndofs_v,nz_qq,nz_vq)

        I_index_qq=Vector{Int}(zeros(nz_qq))
        J_index_qq=Vector{Int}(zeros(nz_qq))
        iptr_elem_qq=Vector{Int}(zeros(Nel))
        S_qq=Vector{Float64}(zeros(nz_qq))
        subset_q=Vector{Int}(zeros(Ndofs_q))

        if Ndofs_v > 0
            I_index_vq=Vector{Int}(zeros(nz_vq))
            J_index_vq=Vector{Int}(zeros(nz_vq))
            iptr_elem_vq=Vector{Int}(zeros(Nel))
            S_vq=Vector{Float64}(zeros(nz_vq))
            S_qvt=Vector{Float64}(zeros(nz_vq))
        else
            I_index_vq=Vector{Int}[]
            J_index_vq=Vector{Int}[]
            iptr_elem_vq=Vector{Int}[]
            S_vq=Vector{Float64}[]
            S_qvt=Vector{Float64}[]
        end


        return new(Ndofs_q, Ndofs_v, I_index_qq, J_index_qq, I_index_vq, J_index_vq,
        iptr_elem_qq, iptr_elem_vq, subset_q, S_qq, S_vq, S_qvt)

    end

end
