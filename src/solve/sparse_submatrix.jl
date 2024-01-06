"""
    sparse_submatrix
"""
function sparse_submatrix(A::SparseMatrixCSC{Int64, Int64}, irows::Vector{Int},icols::Vector{Int})

    irows_out  = Int[]
    icols_out  = Int[]
    vals_out = Float64[]
    II, JJ, Val =findnz(A)
    nvals = size(Val,1)
    for i = 1:nvals
        i_out =  findfirst(x -> x == II[i], irows)
        if i_out != nothing 
            j_out = findfirst(x -> x == JJ[i], icols)
            j_out != nothing  && (irows_out = [irows_out;i_out]; icols_out = [icols_out;j_out]; vals_out = [vals_out; Val[i]])
        end
    end

    ir = size(irows,1) - size(irows_out,1)
    ir > 0   && (println("warning: submatrix has ", ir, " zeros rows eliminated"))
    jr = size(icols,1) - size(icols_out,1)
    jr > 0   && (println("warning: submatrix has ", jr, " zeros columns eliminated"))

    return irows_out,icols_out,vals_out

end 