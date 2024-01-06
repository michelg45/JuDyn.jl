"""
    rigid_body
"""
function rigid_body(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},res::Vector{Float64},
    p::Vector{Float64},alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime::Int,h::Float64)

    S_el = Array{Float64,2}(zeros(12,12))
    res_el  = Vector{Float64}(zeros(12))
    p_el  = Vector{Float64}(zeros(12))

    rbc = Main.SetElements.rigid_body_container
    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container

    iel = findfirst(x -> x == nbr, rbc.numbers)[1]

    inode = rbc.node_orders[iel]
    mass = rbc.masses[iel]
    J = rbc.inertia_tensors[iel]
    
    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed)

    if rotation == true
        t = itime*h
        phi_d = RV3(t*Omega_d)
        grav = rot(-phi_d, mc.gravity)
    else
        grav = mc.gravity
    end
    
    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]
    inv_loc_v = ec.inv_loc_v[iel]

    inv_loc = [inv_loc_x; inv_loc_v .+ mc.Ndofs_q]

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)


    ix = [i for i in 1:3]
    itheta = ix .+ 3
    iv = ix .+ 6
    iomega = itheta .+ 6

    x = cf[inode].x + Vec3(Ds[ix])
    xdot = Vec3(sdot[ix])
    psi = RV3(cf[inode].psi,RV3(Ds[itheta]))
    Dpsi = RV3(Ds[itheta])
    psidot = Vec3(sdot[itheta])
    T = tang(Dpsi)
    W = T*psidot  
    v = Vec3(s[iv]+Ds[iv])
    vdot = Vec3(sdot[iv])
    Omega= Vec3(s[iomega]+Ds[iomega])
    Omegadot= Vec3(sdot[iomega])
    DW = Dtang(psidot,Dpsi)
  
    

    p_el[ix] = mass*grav.v
    


    if rotation == true
        ROmega_d = rot(-psi,Omega_d)
        tOmega_d =  tilde(Omega_d).mat
        res_el[ix] =  mass*(grav - vdot - crossp(Omega_d,v)).v
        res_el[iv]= (v - xdot - crossp(Omega_d,x)).v 
        res_el[itheta] = -(J*Omegadot + crossp(W+ROmega_d, J*Omega)).v
        res_el[iomega]= (Omega - W - ROmega_d).v
        S_el[ix,iv]= mass*(theta_p*eye(3) + alpha_stiff*tOmega_d)
        S_el[iv,ix]= theta_p*eye(3) + alpha_stiff*tOmega_d
        A = (theta_p*I3 + alpha_stiff*tilde(ROmega_d))*T + alpha_stiff*DW
        S_el[itheta,iomega]= (theta_p*J + alpha_stiff*tilde(W+ROmega_d)*J).mat

    else
        res_el[ix] =  mass*(grav-vdot).v
        res_el[iv]= (v - xdot).v
        res_el[itheta] = -(J*Omegadot + crossp(W,J*Omega)).v
        res_el[iomega]= (Omega -W).v
        S_el[ix,iv]= mass*theta_p*eye(3)
        S_el[iv,ix]= theta_p*eye(3)
        A = theta_p*T + alpha_stiff*DW
        S_el[itheta,iomega]= (theta_p*J+alpha_stiff*tilde(W)*J).mat
        
    end

    S_el[iv,iv]=  - alpha_stiff*eye(3)
    S_el[itheta,itheta]= -(tilde(J*Omega)*A).mat
    S_el[iomega,iomega]= - alpha_stiff*eye(3)
    S_el[iomega,itheta]  =  A.mat


    kin_el = 0.5*(mass*dotp(v,v) + dotp(Omega,J*Omega))
    pot_el = - transpose(p_el[1:3])*x.v

    #
    # add element contribution to iteration matrix S and residual vector S
    #
    
    if matrix == false
        push_element_sparse(res,p,iel,inv_loc_x,inv_loc_v,S_el,res_el,p_el)
        return  kin_el, pot_el
    else
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc_x,inv_loc_v,S_el)
    end

end   

"""
    rigid_body_force
"""
function rigid_body_force(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},res::Vector{Float64},p::Vector{Float64},itime::Int,h::Float64)
    
    Ndim = 12
    res_el  = Vector{Float64}(zeros(Ndim))
    p_el  = Vector{Float64}(zeros(Ndim))

    rbc = Main.SetElements.rigid_body_container
    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container

    iel = findfirst(x -> x == nbr, rbc.numbers)[1]

    inode = rbc.node_orders[iel]
    mass = rbc.masses[iel]
    J = rbc.inertia_tensors[iel]
    
    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed)

    if rotation == true
        t = itime*h
        phi_d = RV3(t*Omega_d)
        grav = rot(-phi_d, mc.gravity)
    else
        grav = mc.gravity
    end
    
    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]
    inv_loc_v = ec.inv_loc_v[iel]

    inv_loc = [inv_loc_x; inv_loc_v .+ mc.Ndofs_q]

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)


    ix = [i for i in 1:3]
    itheta = ix .+ 3
    iv = ix .+ 6
    iomega = itheta .+ 6

    x = cf[inode].x + Vec3(Ds[ix])
    xdot = Vec3(sdot[ix])
    psi = RV3(cf[inode].psi,RV3(Ds[itheta]))
    Dpsi = RV3(Ds[itheta])
    psidot = Vec3(sdot[itheta])
    W = tang(Dpsi,psidot)
    v = Vec3(s[iv]+Ds[iv])
    vdot = Vec3(sdot[iv])
    Omega= Vec3(s[iomega]+Ds[iomega])
    Omegadot= Vec3(sdot[iomega])


    p_el[ix] = mass*grav.v
    


    if rotation == true
        ROmega_d = rot(-psi,Omega_d)
        x = cf[inode].x
        res_el[ix] =  mass*(grav - vdot - crossp(Omega_d,v)).v
        res_el[iv]= (v - xdot - crossp(Omega_d,x)).v 
        res_el[itheta] = -(J*Omegadot + crossp(W+ROmega_d, J*Omega)).v
        res_el[iomega]= (Omega - W - ROmega_d).v
    else
        res_el[ix] =  mass*(grav-vdot).v
        res_el[iv]= (v - xdot).v
        res_el[itheta] = -(J*Omegadot + crossp(W,J*Omega)).v
        res_el[iomega]= (Omega -W).v
    end

    kin_el = 0.5*(mass*dotp(v,v) + dotp(Omega,J*Omega))
    pot_el = - transpose(p_el[1:3])*x.v

    #
    # add element contribution to residual vector S
    #
    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    return  kin_el, pot_el

end

