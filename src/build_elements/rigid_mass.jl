"""
    rigid_mass
"""
function rigid_mass(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector,
    alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime,h)

    S_el = Array{Float64,2}(zeros(6,6))
    res_el  = Vector{Float64}(zeros(6))
    p_el  = Vector{Float64}(zeros(6))

    rmc = Main.SetElements.rigid_mass_container
    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container

    iel = findfirst(x -> x == nbr, rmc.numbers)[1]

    inode = rmc.node_orders[iel]
    mass = rmc.masses[iel]
    grav = mc.gravity
    
    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed; tOmega_d = tilde(Omega_d).mat)

    if rotation == true
        phi_d = itime*h*RV3(Omega_d)
        grav = rot(-phi_d, mc.gravity)
    else
        grav = mc.gravity
    end

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x = ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]

    inv_loc = [inv_loc_x; inv_loc_v .+ mc.Ndofs_q]

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

    ix = [i for i in 1:3]
    iv = ix .+ 3

    v = Vec3(s[iv]+Ds[iv])
    vdot = Vec3(sdot[iv])
    xdot = cf[inode].xdot

    Dx = Vec3(Ds[ix])
    x = cf[inode].x  + Dx

    p_el[ix] = mass*grav.v

    rotation == true ? res_el[ix] =  mass*(grav-vdot-crossp(Omega_d,v)).v : res_el[ix] =  mass*(grav-vdot).v
    rotation == true ? res_el[iv]= (v - xdot-crossp(Omega_d,x)).v : res_el[iv]= (v - xdot).v
   
    kin_el = 0.5*mass*dotp(v,v)
    pot_el = - transpose(p_el[1:3])*(cf[inode].x).v

    rotation == true ? S_el[ix,iv]= mass*(theta_p*eye(3)+ alpha_stiff*tOmega_d) : S_el[ix,iv]= mass*theta_p*eye(3)
    rotation == true ? S_el[iv,ix]= theta_p*eye(3)+ alpha_stiff*tOmega_d :  S_el[iv,ix]= theta_p*eye(3)

    S_el[iv,iv]=  - alpha_stiff*eye(3)


    if matrix == false
        push_element_sparse(res,p,iel2,inv_loc_x,inv_loc_v,S_el,res_el,p_el)
        return  kin_el, pot_el
    else
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc_x,inv_loc_v,S_el)
    end
 
end


"""
    rigid_mass_force
"""
function rigid_mass_force(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector,
    itime,h)

    res_el  = Vector{Float64}(zeros(6))
    p_el  = Vector{Float64}(zeros(6))

    rmc = Main.SetElements.rigid_mass_container
    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container

    iel = findfirst(x -> x == nbr, rmc.numbers)[1]

    inode = rmc.node_orders[iel]
    mass = rmc.masses[iel]
    grav = mc.gravity
    
    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed; tOmega_d = tilde(Omega_d).mat)

    if rotation == true
        phi_d = itime*h*RV3(Omega_d)
        grav = rot(-phi_d, mc.gravity)
    else
        grav = mc.gravity
    end

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x = ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]

    inv_loc = [inv_loc_x; inv_loc_v .+ mc.Ndofs_q]

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

    ix = [i for i in 1:3]
    iv = ix .+ 3

    v = Vec3(s[iv]+Ds[iv])
    vdot = Vec3(sdot[iv])
    xdot = cf[inode].xdot
    
    Dx = Vec3(Ds[ix])
    x = cf[inode].x  + Dx

    rotation == true ? res_el[ix] =  mass*(grav-vdot-crossp(Omega_d,v)).v : res_el[ix] =  mass*(grav-vdot).v
    rotation == true ? res_el[iv]= (v - xdot-crossp(Omega_d,x)).v : res_el[iv]= (v - xdot).v
   
    kin_el = 0.5*mass*dotp(v,v)
    pot_el = - transpose(p_el[1:3])*(cf[inode].x).v

    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    return  kin_el, pot_el
    
end