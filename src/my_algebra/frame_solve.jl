
"""
    frame_solve
"""
function frame_solve(H::Vector{NodeFrame},F::Vector,DF::Array,J::Array,Nnodes::Int)
    global tol = 1.e-012
    
    global r = zeros(6)
    global invA = zeros(6,6)
    global invTm = Vector{Matrix}(undef,Nnodes)
    global invT = Vector{Matrix}(undef,Nnodes)
    global p = Vector{Vector}(undef,Nnodes)
   
    ndim = size(DF,2)
    global f = zeros(6,ndim)

    global H_P = NodeFrame(H[1].x,H[1].p)

    niter = 0
    
    nitmax = 20
    test = 1.0e5

    invJ = J\eye(2)

    DF *= invJ

    while (test > tol &&  niter <  nitmax)
   
        r .= 0.0
        invA .= 0.0
        for k = 1:Nnodes
            # p[k] = log_SE3(-(H_0[k]*H_P)*H[k])
            # p[k] = log_SE3(-H_P*H[k]*H_0[k])
            p[k] = log_SE3(-H_P*H[k])
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

    niter == nitmax && (println("max. number of iterations reached in frame_solve")) 

    invA .= 0.0

    

    for k = 1:Nnodes
        invTm[k] = invT_SE3(-p[k])
        invT[k] = invT_SE3(p[k])
        invA += F[k]*invTm[k]
        for nd  = 1:ndim
            f[:,nd] += DF[k,nd]*p[k]
        end
    end

    A = invA\eye(6)
    f = A*f

    Q = Vector{Matrix}(undef,Nnodes)
    for k = 1:Nnodes
        Q[k] = A*(F[k]*invT[k])
    end
    
    g = [rot(H_P.p, Vec3(f[1:3,1])); rot(H_P.p, Vec3(f[1:3,2]))]
    n = crossp(g[1],g[2])/(norm2(g[1])*norm2(g[2]))

    


    global B = Matrix{Matrix}(undef,Nnodes,2)
    for i = 1:2
        for k = 1:Nnodes
            B[k,i] = DF[k,i]*eye(6)+F[k]*DinvT_SE3(-p[k],f[:,i])
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
            E[i] -= B[k,i]*invTm[k]
            R[i,k] = B[k,i]*invT[k]
        end
    end

    # R_0t = eye(6)

       
    for k = 1:Nnodes
        # R_0t[1:3,1:3] = rot(H_0[k].p).mat
        for i = 1:2
            # tr[i] = A*(E[i]*Q[k] + R[i,k])*Adj(-H_0[k])
            tr[i] = A*(E[i]*Q[k] + R[i,k])
            # tr[i] = A*(E[i]*Q[k] + R[i,k])*R_0t
        end
        k == 1 ? QN = vcat(tr[1], tr[2]) : QN = hcat(QN, vcat(tr[1], tr[2]))
    end
    return f, QN
end


