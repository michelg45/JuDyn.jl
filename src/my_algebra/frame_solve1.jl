"""
    frame_solve1
"""
function frame_solve1(H::Vector{NodeFrame},F::Vector,DF::Array,Nnodes::Int)

    global p = Vector{Vector}(undef,Nnodes)

    ndim = size(DF,2)
    global f_0 = zeros(6,ndim)
    global r = zeros(6)
    global invA = zeros(6,6)

    global tol = 1.e-012
    global H_P = H[1]

    
   
    Nitmax = 20
    

    niter = 0    
    test = 1.0e5

    while (test > tol &&  niter <  Nitmax)
   
        r .= 0.0
        invA .= 0.0
        for k = 1:Nnodes
            # p[k] = log_SE3(-(H_P)*H[k]*H_0[k])
            p[k] = log_SE3(-H_P*H[k])
            # p[k] = log_SE3(-(H_0[k]*H_P)*H[k])
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

    invA .= 0.0
    
    for k = 1:Nnodes
        invA += F[k]*invT_SE3(-p[k])
        for nd  = 1:ndim
            f_0[:,nd] += DF[k,nd]*p[k]
        end
    end
    f_0 = invA\f_0

    J = f_0[1:2,1:2]

    invJ = J\eye(2)

    

    g_0 = [rot(H_P.p, Vec3(f_0[1:3,1])); rot(H_P.p, Vec3(f_0[1:3,2]))]
    n_0 = crossp(g_0[1],g_0[2])/(norm2(g_0[1])*norm2(g_0[2]))

    f_0 = f_0*invJ


    return H_P, f_0, g_0, n_0, J

end