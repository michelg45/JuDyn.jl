function beam(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64}, alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime::Int,h::Float64)

    bc = SetElements.beam_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    mc =  Main.model_container

    iel = findfirst(x -> x == nbr, bc.numbers)[1]

    constant_inertia = bc.constant_inertia[iel]

    Ndim = 24-Int(constant_inertia)*6

    S_el = Array{Float64,2}(zeros(Ndim,Ndim))
    res_el  = Vector{Float64}(zeros(Ndim))

    
    DstrDq = Array{Float64,2}(zeros(6,12))
    DstrDp = Array{Float64,2}(zeros(6,12))

    Ndofs_q = mc.Ndofs_q

    

    Z33 = zeros(3,3)

    inodes = bc.node_orders[iel]
    Length = bc.length[iel]
    psi_rel = bc.local_node_orientations[iel]
    R_rel = [rot(psi_rel[i]) for i in 1:2]
    stiffness_properties = bc.stiffness_properties[iel]
    mass_properties = bc.mass_properties[iel]
    mass = mass_properties[1]*Length

    iner_rot = Length*diagm(mass_properties[2:4])
    iner_tr = mass*eye(3)

    if constant_inertia == false
        mass_kernel_tr = 1.0/6.0*iner_tr
        Gyr_kernel = mass_properties[1]*[1.0/3.0 1.0/6.0; 1.0/6.0 1.0/3.0]
        mass_kernel_rot =  1.0/6.0*(iner_rot)
        mass_tr = [2*mass_kernel_tr mass_kernel_tr; mass_kernel_tr 2*mass_kernel_tr]
        mass_rot = [2*mass_kernel_rot mass_kernel_rot; mass_kernel_rot 2*mass_kernel_rot]
    end

    Trot = [R_rel[1].mat Z33; Z33 R_rel[2].mat]
    
    constant_inertia == true ? (U6 = [eye(3); eye(3)]; mass_rot = 0.5* transpose(Trot)*U6*iner_rot) : mass_rot = transpose(Trot)*mass_rot*Trot
        
    rotation = mc.uniform_rotation

    rotation == true && (Omega_d = mc.rotation_speed; tOmega_d = tilde(Omega_d).mat; Omega_d = Omega_d.v;
    Gyr_x = [tOmega_d zeros(3,3); zeros(3,3) tOmega_d])

    if rotation == true
        t = itime*h
        phi_d = RV3(t*Omega_d)
        gravity = (rot(-phi_d, mc.gravity)).v
    else
        gravity = (mc.gravity).v 
    end
    # 

    x  =  Vector{Vec3}(undef,2)
    Dx  =  Vector{Vec3}(undef,2)
    psi  =  Vector{RV3}(undef,2)
    Dpsi  =  Vector{RV3}(undef,2)
    T  = Vector{Mat3}(undef,2)
    W  = Vector{Vec3}(undef,2)
    DW  = Vector{Mat3}(undef,2)



    for i = 1:2
        x[i] = cf[inodes[i]].x
        Dx[i] = cf[inodes[i]].Dx
        psi[i] = cf[inodes[i]].psi
        Dpsi[i] = cf[inodes[i]].Dpsi
        T[i] = cf[inodes[i]].T
        W[i] = cf[inodes[i]].W
        DW[i] = cf[inodes[i]].DW
    end


    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]
    inv_loc= [inv_loc_x; inv_loc_v.+Ndofs_q]




    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)


    R_1 = rot(psi[1])
    phi = RV3(-RV3(psi[1],-psi_rel[1]),RV3(psi[2],-psi_rel[2]))
    u = rot(-psi[1], x[2]-x[1])
    tdu = tilde(u)
    u = rot(psi_rel[1],u)
    invT = invtang(phi)
    tinvTR = transpose(invT)*R_rel[1]
    DinvTTu = DinvtangT(u,phi)
    invTRT = 1.0/Length*(tinvTR*transpose(R_1)).mat
    A2 = 1.0/Length*DinvTTu*invT*R_rel[2]
    A1 = 1.0/Length*(tinvTR*tdu - DinvTTu*tinvTR)

    e1 = zeros(3)
    e1[1] = 1.0
    strain = [1/Length*(transpose(invT)*u).v - e1; 1.0/Length*phi.v]
    stress = stiffness_properties .* strain

    DstrDp =    [-invTRT  A1.mat  invTRT  A2.mat; Z33 -1.0/Length*(tinvTR).mat Z33 1.0/Length*(invT*R_rel[2]).mat]
    DstrDq =    [-invTRT  (A1*T[1]).mat  invTRT  (A2*T[2]).mat; Z33 -1.0/Length*(tinvTR*T[1]).mat Z33 1.0/Length*(invT*R_rel[2]*T[2]).mat]


    idepl = [i for i in 1:12]
    ix = [[i for i in 1:3]; [i for i in 7:9]]
    itheta = ix .+ 3
    constant_inertia == true ?  ivx = [i for i in 13:15] : ivx= ix .+ 12

    iomega = ivx .+3

    xdot = sdot[ix]
    vdot = sdot[ivx]
    v = s[ivx] + Ds[ivx]
    Omega= s[iomega]+Ds[iomega]
    Omegadot= sdot[iomega]

    tW = [tilde(W[1]).mat Z33; Z33 tilde(W[2]).mat]
    vW = [W[1].v; W[2].v]

    res_el[idepl] = - Length*transpose(DstrDp)*stress

    S_el[idepl,idepl] = alpha_stiff*Length*transpose(DstrDp)*diagm(stiffness_properties)*DstrDq 
    
    p_el = 0.5*mass_properties[1]*Length*[gravity;  gravity]

    if rotation == true 
        t2Omega_d = [tOmega_d Z33; Z33 tOmega_d]
        x_tr = [x[1].v; x[2].v]
        f_iner_tr =  mass_tr*(vdot+ Gyr_x*v) 
        f_iner_rot = mass_rot*Omegadot + (tW + t2Omega_d)*mass_rot*Omega
        S_el[ix,ivx] = mass_tr*(theta_p*eye(6) + alpha_stiff*Gyr_x)
        S_el[itheta,iomega] = (theta_p*eye(6) + alpha_stiff*(t2Omega_d + tW))*mass_rot
        res_el[ivx] =  v - xdot - Gyr_x*x_tr
        S_el[ivx,ix] = theta_p*eye(6) +  alpha_stiff*Gyr_x
        res_el[iomega] =  Omega - vW - [Omega_d; Omega_d]
    else
        
        f_iner_rot = mass_rot*Omegadot + tW*mass_rot*Omega
        res_el[itheta] -= f_iner_rot
        S_el[itheta,iomega] = (theta_p*eye(6) + alpha_stiff*tW)*mass_rot
        
    end

    global A = [(theta_p*T[1]+alpha_stiff*DW[1]).mat  Z33; Z33 (theta_p*T[2]+alpha_stiff*DW[2]).mat]

        if constant_inertia == true 
            TrW = Trot*vW
            f_iner_tr = 0.5*mass*[vdot; vdot]
            S_el[ix,ivx] = 0.5*mass*theta_p*U6
            res_el[ivx] =  v - 0.5*(xdot[1:3]+xdot[4:6])
            S_el[ivx,ix] = 0.5*theta_p*U6'
            res_el[iomega] =  Omega - 0.5*(TrW[1:3] + TrW[4:6])
            S_el[ivx,ivx] = -alpha_stiff*eye(3)
            S_el[iomega,iomega] = -alpha_stiff*eye(3)
            S_el[iomega, itheta] = 0.5*(A[1:3,1:6]+A[4:6,1:6])
        else
            S_el[ix,ivx] = theta_p*mass_tr
            f_iner_tr = mass_tr*vdot 
            res_el[ivx] =  v - xdot
            S_el[ivx,ix] = theta_p*eye(6)
            res_el[iomega] =  Omega - vW
            S_el[ivx,ivx] = -alpha_stiff*eye(6)
            S_el[iomega,iomega] = -alpha_stiff*eye(6)
            S_el[iomega, itheta] = A
        end
    

    
    res_el[ix] +=   p_el - f_iner_tr


    str_el = 0.5*Length*transpose(strain)*stress
    constant_inertia == true ? kin_el = 0.5*mass*(transpose(v)*v +  transpose(Omega)*iner_rot*Omega) : kin_el = 0.5*mass*(transpose(v)*v +  transpose(Omega)*mass_rot*Omega)

    if matrix == false
        if norm(gravity) == 0.0
            pot_el = 0.0
            push_element_sparse(res,iel2,inv_loc_x,inv_loc_v,S_el,res_el)
        else
            pot_el = p_el'*[x[1].v; x[2].v]
            p_el = [p_el[1:3]; zeros(3); p_el[4:6]; zeros(3)]
            push_element_sparse(res,p,iel2,inv_loc_x,inv_loc_v,S_el,res_el,p_el)
        end

        return kin_el, str_el, pot_el
        
    else
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc_x,inv_loc_v,S_el)
    end

end

function beam(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64},theta_p::Float64,itime::Int,h::Float64)
    alpha_stiff = 1.0
    matrix = false

    return kin_el, str_el, pot_el = beam(nbr,y,Dy,ydot,
        res,p, alpha_stiff,theta_p,matrix,itime::Int,h::Float64)
end

function beam(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},matrix_type::String,itime::Int,h::Float64)
    
    matrix_type == "stiffness" && (alpha_stiff = 1.0; theta_p = 0.0)
    matrix_type == "mass" && (alpha_stiff = 0.0; theta_p = 1.0)
    matrix = true
    
    N = size(y,1)
    res = zeros(N)
    p = zeros(N)


    return inv_loc, S_el = beam(nbr,y,Dy,ydot, res,p, alpha_stiff,theta_p,matrix,itime,h)
end
