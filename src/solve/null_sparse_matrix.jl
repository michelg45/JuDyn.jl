"""
    null_sparse_matrix
"""
function null_sparse_matrix(::SparseMatrixStructure)
    sps = Main.Solve.sparse_matrix_S
    sps.S_qq[:] .= 0.0
    sps.S_vq[:] .= 0.0
    sps.S_qvt[:] .= 0.0 

end
