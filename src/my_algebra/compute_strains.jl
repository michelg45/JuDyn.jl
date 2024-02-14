"""
    compute_strains
"""
function compute_strains(H::Vector{NodeFrame},H_0::Vector{NodeFrame},p::Vector{Vector},F::Vector,DF::Array,J::Array,Nnodes::Int)

    global invA = zeros(6,6)
    global invTm = Vector{Matrix}(undef,Nnodes)
    global invT = Vector{Matrix}(undef,Nnodes)
      
    ndim = size(DF,2)
    global f = zeros(6,ndim)

    invJ = J\eye(2)


    # println("frame_solve2 niter ", niter)
    # niter == nitmax && (println("max. number of iterations reached in frame_solve2")) 

    invA .= 0.0
    
    for k = 1:Nnodes
        invTm[k] = invT_SE3(-p[k])
        invT[k] = invT_SE3(p[k])
        invA += F[k]*invTm[k]
        for nd  = 1:ndim
            f[:,nd] += DF[k,nd]*p[k]
        end
    end
    f = (invA\f)*invJ

    B = Vector{Matrix}(undef,Nnodes)
    for k = 1:Nnodes
        B[k] = invA\(F[k]*invT[k]*Adj(-H_0[k]))
    end
    
    DF *= invJ


    global C = Matrix{Matrix}(undef,Nnodes,2)
    for i = 1:2
        for k = 1:Nnodes
            C[k,i] = DF[k,i]*eye(6)+F[k]*DinvT_SE3(-p[k],f[:,i])
        end
    end

    global E = Vector{Matrix}(undef,2)
    global R = Matrix{Matrix}(undef, 2, Nnodes)
    global tr = Matrix{Matrix}(undef, 2, 1)
    global QN
    global Q


    for i = 1:2
        E[i] = zeros(6,6)
        for k = 1:Nnodes
            E[i] -= C[k,i]*invTm[k]
            R[i,k] = C[k,i]*invT[k]*Adj(-H_0[k])
        end
    end
    
    for k = 1:Nnodes

        k == 1 ? Q = B[1] : Q = hcat(Q,B[k])

        rot_k = rot(H[k].p).mat
        for i = 1:2
            tr[i] = invA\(E[i]*B[k] + R[i,k])
        #    tr[i][1:3,1:3] = tr[i][1:3,1:3]*rot_k
        end
        k == 1 ? QN = vcat(tr[1], tr[2]) : QN = hcat(QN, vcat(tr[1], tr[2]))
    end

    return f, QN
end


