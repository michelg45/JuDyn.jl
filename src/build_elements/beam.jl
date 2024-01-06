"""
    beam

    function constructing the residual vector and the iteration matrix of a beam element.

    Depending on the visco_type::String poarameter (provided by beam_container), the material can be elastic or viscoelastic (Maxwell or Kevin-Voigt).
        * visco_type = "none"       elastic material
        * visco_type = "maxwell"    2-branch Maxwell viscoelastic material
        * visco_type = "damped"     1-branch Kevin-Voigt viscoelastic material
        
        calling sequence:

        kin_el, str_el, pot_el = beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)


"""
function beam(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64}, alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime::Int,h::Float64,niter::Int64)

    bc = SetElements.beam_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    mc =  Main.model_container

    iel = findfirst(x -> x == nbr, bc.numbers)[1]

    S_el = Array{Float64,2}(zeros(24,24))
    res_el  = Vector{Float64}(zeros(24))

    
    DstrDq = Array{Float64,2}(zeros(6,12))
    DstrDp = Array{Float64,2}(zeros(6,12))

    Ndofs_q = mc.Ndofs_q

    Z33 = zeros(3,3)

    inodes = bc.node_orders[iel]
    Length = bc.length[iel]
    psi_rel = bc.local_node_orientations[iel]
    R_rel = [rot(psi_rel[i]) for i in 1:2]

    K = bc.stiffness_properties[iel]
    mass_properties = bc.mass_properties[iel]
    mass = mass_properties[1]*Length

    visco_type = bc.visco_type[iel]

    if visco_type == "maxwell"
    
        time_constants = bc.time_constants[iel]
        ratio_infty = bc.ratio_infty[iel]
        K_b = (1.0 - ratio_infty)*K

        if itime > 1
            Gamma_2 = broadcast(x -> h/x/2.0*exp(-h/2.0/x), time_constants)
            Gamma_1 = broadcast(x -> exp(-h/x), time_constants)
        else
            Gamma_1 = zeros(6)
            Gamma_2 = zeros(6)
        end

        niter == 1 && (bc.strains[iel][:,1] = bc.strains[iel][:,2];
        bc.visco_strains[iel][:,1] = bc.visco_strains[iel][:,2])

end


    iner_rot = Length*diagm(mass_properties[2:4])
    iner_tr = mass*eye(3)

    mass_kernel_tr = 1.0/6.0*iner_tr
    Gyr_kernel = mass_properties[1]*[1.0/3.0 1.0/6.0; 1.0/6.0 1.0/3.0]
    mass_kernel_rot =  1.0/6.0*(iner_rot)
    mass_tr = [2*mass_kernel_tr mass_kernel_tr; mass_kernel_tr 2*mass_kernel_tr]
    mass_rot = [2*mass_kernel_rot mass_kernel_rot; mass_kernel_rot 2*mass_kernel_rot]
    
    Trot = [R_rel[1].mat Z33; Z33 R_rel[2].mat]
    mass_rot = transpose(Trot)*mass_rot*Trot
        
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

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]
    inv_loc= [inv_loc_x; inv_loc_v.+Ndofs_q]

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

    idepl = [i for i in 1:12]
    ix = [[i for i in 1:3]; [i for i in 7:9]]
    itheta = ix .+ 3
    ivx= ix .+ 12

    iomega = ivx .+3

    xdot = sdot[ix]
    vdot = sdot[ivx]
    v = s[ivx] + Ds[ivx]
    Omega= s[iomega]+Ds[iomega]
    Omegadot= sdot[iomega]

    Dx_1 = Vec3(Ds[1:3])
    psidot_1 = Vec3(sdot[4:6])
    Dx_2 = Vec3(Ds[7:9])
    psidot_2 = Vec3(sdot[10:12])
    Dpsi_1 = RV3(Ds[4:6])
    Dpsi_2 = RV3(Ds[10:12])
    x_1 = cf[inodes[1]].x + Dx_1
    x_2 = cf[inodes[2]].x + Dx_2
    
    psi_1 = RV3(cf[inodes[1]].psi, Dpsi_1) 
    psi_2 = RV3(cf[inodes[2]].psi, Dpsi_2)
    R_1 = rot(psi_1) 

    T_1 = tang(Dpsi_1)
    W_1 = tang(Dpsi_1)*psidot_1
    DW_1 = Dtang(psidot_1,Dpsi_1)
    T_2 = tang(Dpsi_2)
    W_2 = tang(Dpsi_2)*psidot_2
    DW_2 = Dtang(psidot_2,Dpsi_2)
   
    # 
    phi = RV3(-RV3(psi_1,-psi_rel[1]),RV3(psi_2,-psi_rel[2]))
    u = rot(-psi_1, x_2-x_1)
    tdu = tilde(u)
    u = rot(-psi_rel[1],u)
    invT = invtang(phi)


"""    tinvTR = transpose(invT)*R_rel[1]
    DinvTTu = DinvtangT(u,phi)
    invTRT = 1.0/Length*(tinvTR*transpose(R_1)).mat
    A2 = 1.0/Length*DinvTTu*invT*R_rel[2]
    A1 = 1.0/Length*(tinvTR*tdu - DinvTTu*tinvTR)
"""
    tinvTR = transpose(invT)*transpose(R_rel[1])
    DinvTTu = DinvtangT(u,phi)
    invTRT = 1.0/Length*(tinvTR*transpose(R_1)).mat
    A2 = 1.0/Length*DinvTTu*invT*transpose(R_rel[2])
    A1 = 1.0/Length*(tinvTR*tdu - DinvTTu*tinvTR)

    e1 = zeros(3)
    e1[1] = 1.0
    strain = [1/Length*(transpose(invT)*u).v - e1; 1.0/Length*phi.v]

    if visco_type == "maxwell"
        alpha = Gamma_1 .* bc.visco_strains[iel][:,1] + Gamma_2 .* (strain + bc.strains[iel][:,1]) 
        stress =  K .* strain - K_b .* alpha
        bc.visco_strains[iel][:,2] = alpha
    else
        stress = K .* strain
    end

    bc.stresses[iel][:] = stress
    bc.strains[iel][:,2] = strain

    """DstrDp =    [-invTRT  A1.mat  invTRT  A2.mat; Z33 -1.0/Length*(tinvTR).mat Z33 1.0/Length*(invT*R_rel[2]).mat]
    DstrDq =    [-invTRT  (A1*T_1).mat  invTRT  (A2*T_2).mat; Z33 -1.0/Length*(tinvTR*T_1).mat Z33 1.0/Length*(invT*R_rel[2]*T_2).mat]
"""
    DstrDp =    [-invTRT  A1.mat  invTRT  A2.mat; Z33 -1.0/Length*(tinvTR).mat Z33 1.0/Length*(invT*transpose(R_rel[2])).mat]
    DstrDq =    [-invTRT  (A1*T_1).mat  invTRT  (A2*T_2).mat; Z33 -1.0/Length*(tinvTR*T_1).mat Z33 1.0/Length*(invT*transpose(R_rel[2])*T_2).mat]




    tW = [tilde(W_1).mat Z33; Z33 tilde(W_2).mat]
    vW = [W_1.v; W_2.v]

    res_el[idepl] = - Length*transpose(DstrDp)*stress



    S_el[idepl,idepl] = alpha_stiff*Length*transpose(DstrDp)*diagm(K)*DstrDq 
    
    if visco_type == "damped"
        damping_properties = bc.time_constants[iel] .* K
        qdot = sdot[idepl]
        res_el[idepl] -= Length*transpose(DstrDp)*diagm(damping_properties)*DstrDq*qdot
        S_el[idepl,idepl] += theta_p*Length*transpose(DstrDp)*diagm(damping_properties)*DstrDq 
    end 

    p_el = 0.5*mass_properties[1]*Length*[gravity;  gravity ]

    if rotation == true 
        t2Omega_d = [tOmega_d Z33; Z33 tOmega_d]
        x_tr = [x_1.v; x_2.v]
        f_iner_tr =  mass_tr*(vdot+ Gyr_x*v) 
        f_iner_rot = mass_rot*Omegadot + (tW + t2Omega_d)*mass_rot*Omega
        S_el[ix,ivx] = mass_tr*(theta_p*eye(6) + alpha_stiff*Gyr_x)
        S_el[itheta,iomega] = (theta_p*eye(6) + alpha_stiff*(t2Omega_d + tW))*mass_rot
        res_el[ivx] =  v - xdot - Gyr_x*x_tr
        S_el[ivx,ix] = theta_p*eye(6) +  alpha_stiff*Gyr_x
        res_el[iomega] =  Omega - vW - [Omega_d; Omega_d]
    else

        f_iner_tr = mass_tr*vdot
        f_iner_rot = mass_rot*Omegadot + tW*mass_rot*Omega
        S_el[ix,ivx] = theta_p*mass_tr
        S_el[itheta,iomega] = (theta_p*eye(6) + tW)*mass_rot
        res_el[ivx] =  v - xdot
        S_el[ivx,ix] = theta_p*eye(6)
        res_el[iomega] =  Omega - vW
    end       
        
    res_el[itheta] -= f_iner_rot
    S_el[itheta,iomega] = (theta_p*eye(6) + alpha_stiff*tW)*mass_rot
        
    

    A = [(theta_p*T_1+alpha_stiff*DW_1).mat  Z33; Z33 (theta_p*T_2+alpha_stiff*DW_2).mat]


    S_el[ivx,ivx] = -alpha_stiff*eye(6)
    S_el[iomega,iomega] = -alpha_stiff*eye(6)
    S_el[iomega, itheta] = [(theta_p*T_1+DW_1).mat  Z33; Z33 (theta_p*T_2+DW_2).mat]
    res_el[ix] +=   p_el - f_iner_tr


 

    str_el = 0.5*Length*transpose(strain)*stress
    kin_el = 0.5*(transpose(v)*mass_tr*v +  transpose(Omega)*mass_rot*Omega)

    if matrix == false
        if norm(gravity) == 0.0
            pot_el = 0.0
            push_element_sparse(res,iel2,inv_loc_x,inv_loc_v,S_el,res_el)
        else
            pot_el = p_el'*[x_1.v; x_2.v]
            p_el = [p_el[1:3]; zeros(3); p_el[4:6]; zeros(3)]
            push_element_sparse(res,p,iel2,inv_loc_x,inv_loc_v,S_el,res_el,p_el)
        end

        return kin_el, str_el, pot_el
        
    else
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc_x,inv_loc_v,S_el)
    end

end


"""
    beam_force

    function constructing the residual vector  of a beam element.

        Depending on the visco_type::String poarameter (provided by beam_container), the material can be elastic or viscoelastic (Maxwell or Kevin-Voigt).
            * visco_type = "none"       elastic material
            * visco_type = "maxwell"    2-branch Maxwell viscoelastic material
            * visco_type = "damped"     1-branch Kevin-Voigt viscoelastic material
        
        calling sequence:

        kin_el, str_el, pot_el = beam_force(nbr,y_n,Dy,ydot_np1,res,p,itime,h)

        
"""
function beam_force(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64},itime::Int,h::Float64,niter::Int64)

    bc = SetElements.beam_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    mc =  Main.model_container

    iel = findfirst(x -> x == nbr, bc.numbers)[1]

    Ndim = 24
    res_el  = Vector{Float64}(zeros(Ndim))

    DstrDp = Array{Float64,2}(zeros(6,12))

    Ndofs_q = mc.Ndofs_q

    visco_type = bc.visco_type[iel]

    Z33 = zeros(3,3)

    inodes = bc.node_orders[iel]
    Length = bc.length[iel]
    psi_rel = bc.local_node_orientations[iel]
    R_rel = [rot(psi_rel[i]) for i in 1:2]
    K = bc.stiffness_properties[iel]
    mass_properties = bc.mass_properties[iel]
    mass = mass_properties[1]*Length

    if visco_type == "maxwell"
    
        time_constants = bc.time_constants[iel]
        ratio_infty = bc.ratio_infty[iel]
        K_b = (1.0 - ratio_infty)*K

        if itime > 1
            Gamma_2 = broadcast(x -> h/x/2.0*exp(-h/2.0/x), time_constants)
            Gamma_1 = broadcast(x -> exp(-h/x), time_constants)
        else
            Gamma_1 = zeros(6)
            Gamma_2 = zeros(6)
        end

        niter == 1 && (bc.strains[iel][:,1] = bc.strains[iel][:,2];
        bc.visco_strains[iel][:,1] = bc.visco_strains[iel][:,2])

end

    iner_rot = Length*diagm(mass_properties[2:4])
    iner_tr = mass*eye(3)

    mass_kernel_tr = 1.0/6.0*iner_tr
    mass_kernel_rot =  1.0/6.0*(iner_rot)
    mass_tr = [2*mass_kernel_tr mass_kernel_tr; mass_kernel_tr 2*mass_kernel_tr]
    mass_rot = [2*mass_kernel_rot mass_kernel_rot; mass_kernel_rot 2*mass_kernel_rot]
    
    Trot = [R_rel[1].mat Z33; Z33 R_rel[2].mat]
    mass_rot = transpose(Trot)*mass_rot*Trot
        
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

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]
    inv_loc= [inv_loc_x; inv_loc_v.+Ndofs_q]


    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

    idepl = [i for i in 1:12]
    ix = [[i for i in 1:3]; [i for i in 7:9]]
    itheta = ix .+ 3
    ivx= ix .+ 12

    iomega = ivx .+3

    xdot = sdot[ix]
    vdot = sdot[ivx]
    v = s[ivx] + Ds[ivx]
    Omega= s[iomega]+Ds[iomega]
    Omegadot= sdot[iomega]

    Dx_1 = Vec3(Ds[1:3])
    psidot_1 = Vec3(sdot[4:6])
    Dx_2 = Vec3(Ds[7:9])
    psidot_2 = Vec3(sdot[10:12])
    Dpsi_1 = RV3(Ds[4:6])
    Dpsi_2 = RV3(Ds[10:12])
    x_1 = cf[inodes[1]].x + Dx_1
    x_2 = cf[inodes[2]].x + Dx_2

    psi_1 = RV3(cf[inodes[1]].psi, Dpsi_1) 
    psi_2 = RV3(cf[inodes[2]].psi, Dpsi_2)
    R_1 = rot(psi_1) 

    W_1 = tang(Dpsi_1)*psidot_1
    W_2 = tang(Dpsi_2)*psidot_2

    
    phi = RV3(-RV3(psi_1,-psi_rel[1]),RV3(psi_2,-psi_rel[2]))
    u = rot(-psi_1, x_2-x_1)
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

    if visco_type == "maxwell"
        alpha = Gamma_1 .* bc.visco_strains[iel][:,1] + Gamma_2 .* (strain + bc.strains[iel][:,1]) 
        stress =  K .* strain - K_b .* alpha
        bc.visco_strains[iel][:,2] = alpha
    else
        stress = K .* strain
    end

    bc.stresses[iel][:] = stress
    bc.strains[iel][:,2] = strain

    DstrDp =    [-invTRT  A1.mat  invTRT  A2.mat; Z33 -1.0/Length*(tinvTR).mat Z33 1.0/Length*(invT*R_rel[2]).mat]
    

    tW = [tilde(W_1).mat Z33; Z33 tilde(W_2).mat]
    vW = [W_1.v; W_2.v]

    res_el[idepl] = - Length*transpose(DstrDp)*stress

    if bc.visco_type[iel] == "damped"
        damping_properties = bc.time_constants[iel] .* K
        qdot = sdot[idepl]
        T_1 = tang(Dpsi_1)
        T_2 = tang(Dpsi_2)
        DstrDq =    [-invTRT  (A1*T_1).mat  invTRT  (A2*T_2).mat; Z33 -1.0/Length*(tinvTR*T_1).mat Z33 1.0/Length*(invT*R_rel[2]*T_2).mat]
        res_el[idepl] -= Length*transpose(DstrDp)*diagm(damping_properties)*DstrDq*qdot
    end

      
    p_el = 0.5*mass_properties[1]*Length*[gravity; gravity]

    if rotation == true 
        t2Omega_d = [tOmega_d Z33; Z33 tOmega_d]
        x_tr = [x_1.v; x_2.v]
        f_iner_tr =  mass_tr*(vdot+ Gyr_x*v) 
        f_iner_rot = mass_rot*Omegadot + (tW + t2Omega_d)*mass_rot*Omega
        res_el[ivx] =  v - xdot - Gyr_x*x_tr
        res_el[iomega] =  Omega - vW - [Omega_d; Omega_d]
    else

        f_iner_tr = mass_tr*vdot
        f_iner_rot = mass_rot*Omegadot + tW*mass_rot*Omega
        res_el[ivx] =  v - xdot
        res_el[iomega] =  Omega - vW
    end       
        
    res_el[itheta] -= f_iner_rot
    res_el[ix] +=   p_el - f_iner_tr

    str_el = 0.5*Length*transpose(strain)*stress
    kin_el = 0.5*(transpose(v)*mass_tr*v +  transpose(Omega)*mass_rot*Omega)

    pot_el = p_el'*[x_1.v; x_2.v]
    p_el = [p_el[1:3]; zeros(3); p_el[4:6]; zeros(15)]

    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i];  p[iloc] += p_el[i]) 
    end



    return kin_el, str_el, pot_el
        

end

