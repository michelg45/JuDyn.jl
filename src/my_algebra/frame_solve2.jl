function frame_solve2(H::Vector{NodeFrame},H_0::Vector{NodeFrame},F::Vector,DF::Array,J::Array,Nnodes::Int)
    global tol = 1.e-09
    
    global r = zeros(6)
    global invA = zeros(6,6)
    global invTm = Vector{Matrix}(undef,Nnodes)
    global invT = Vector{Matrix}(undef,Nnodes)
    global p = Vector{Vector}(undef,Nnodes)
   
    ndim = size(DF,2)
    global f = zeros(6,ndim)

    global H_P = NodeFrame(H[1].x,H[1].p)

    niter = 0
    
    nitmax = 10
    test = 1.0e5

    invJ = J\eye(2)

    while (test > tol &&  niter <  nitmax)
   
        r .= 0.0
        invA .= 0.0
        for k = 1:Nnodes
            p[k] = log_SE3(-H_P*NodeFrame(H[k].x, RV3(H[k].p, H_0[k].p)))
            r += F[k]*p[k]
        end

        test = norm(r)
        niter  == 0 && (tol = tol/test)
        for k = 1:Nnodes
            invA += F[k]*invT_SE3(-p[k])
        end
        dh = invA\r
        H_P = H_P*exp_SE3(dh)
        niter += 1
    end
    # println("frame_solve2 nit_solve ", niter)
    # println("p ", p)
    # println("H_P ", H_P)
    niter == nitmax && (println("max. number of iterations reached in frame_solve2")) 

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
    
    g = [rot(H_P.p, Vec3(f[1:3,1])); rot(H_P.p, Vec3(f[1:3,2]))]
    n = crossp(g[1],g[2])/(norm2(g[1])*norm2(g[2]))

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

    # println("f ", f)

    return f, QN
end


