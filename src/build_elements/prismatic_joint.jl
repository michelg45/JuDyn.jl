"""
    prismatic_joint
"""
function prismatic_joint(nbr::Int,Dy::Vector,y_n::Vector,ydot::Vector,res::Vector,p::Vector,alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime::Int,h::Float64)


    pjc = SetElements.prism_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, pjc.number)[1]

    str_time_function = pjc.time_function[iel]

    params = pjc.params[iel]
    k = pjc.scale_factor[iel]
    RV_rel = pjc.RV_rel[iel]
    RV_line = pjc.RV_line[iel]
    l_0 =  pjc.l_0[iel]
    
    mode = pjc.mode[iel]
    nodes =  pjc.node_orders[iel]

    axis = rot(RV_line,Vec3(1.0,0.0,0.0))
    

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    mode == "driven" ? ndim = 18 : ndim =19


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    nodes = pjc.node_orders[iel]
    inode1 = nodes[1]
    inode2 = nodes[2]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult= ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]

    inv_loc = [inv_loc_x; inv_loc_s]

    str_time_function != " " &&  (disp = input_functions[str_time_function])

    s, Ds, sdot = pull_vectors(inv_loc, y_n, Dy, ydot)

    i_x = [i for i=1:12]
    i_x1 = i_x[1:3]
    i_x2 = i_x[7:9]
    itheta_1 = i_x[4:6]
    itheta_2 = i_x[10:12]

    if mode != "driven"
        i_int = 13
        imult= [i for i=14:19]
        x_int_n = s[i_int]
        xdot_int = sdot[i_int]
        Dx_int = Ds[i_int]
    else
        imult = [i for i=13:18]
    end

    imult_x = imult[1:3]
    imult_theta = imult[4:6]
    mult_x = Vec3(Ds[imult_x])
    mult_x_n = s[imult_x]
    mult_x_av = 0.5*(mult_x.v + mult_x_n)
    mult_theta = Vec3(Ds[imult_theta])

    Dpsi_1 = RV3(Ds[itheta_1])
    Dpsi_2 = RV3(Ds[itheta_2])
    x_1 = cf[inode1].x + Vec3(Ds[i_x1])
    x_2 = cf[inode2].x + Vec3(Ds[i_x2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)

    T_1 = (tang(Dpsi_1)).mat
    T_2 = (tang(Dpsi_2)).mat
    R_1 = (rot(psi_1)).mat

    u_psi = RV3(-RV_rel,RV3(-psi_1, psi_2))

    if mode == "driven"
        elong  = time_function(itime,h,params)
        delong  = elong - time_function(itime-1,h,params)
    else
        elong = x_int_n + Dx_int
        delong = Dx_int
    end

    u_psi = RV3(-RV_rel,RV3(-psi_1, psi_2))
    X = (l_0+elong)*axis  
    u_x = x_2 - x_1 - rot(psi_1,X)

    rtmult_x = rot(-psi_1,mult_x)

    res_el[imult_x] = k*u_x.v
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x1] = -k*mult_x.v
    res_el[i_x2] = k*mult_x.v
    res_el[itheta_1] =  - k*rot(RV_rel,mult_theta).v - k*(crossp(X,rtmult_x)).v
    res_el[itheta_2] =    k*mult_theta.v

    R_rel = rot(RV_rel).mat

    S_el[i_x1, imult_x] = alpha_stiff*k*eye(3)
    S_el[i_x2, imult_x] = - alpha_stiff*k*eye(3)
    S_el[imult_x,i_x1] = alpha_stiff*k*(eye(3))
    S_el[imult_x,i_x2] = - alpha_stiff*k*(eye(3))
    S_el[imult_x,itheta_1] = -alpha_stiff*k*R_1*tilde(X).mat*T_1
    S_el[itheta_1,itheta_1] =   alpha_stiff*k*((tilde(X)*tilde(rtmult_x)).mat)*T_1
    S_el[itheta_1,imult_x] = alpha_stiff*k*tilde(X).mat*R_1'
    
    S_el[itheta_1,imult_theta] =   alpha_stiff*k*R_rel
    S_el[itheta_2,imult_theta] =   -alpha_stiff*k*eye(3)
    S_el[imult_theta,itheta_1] = alpha_stiff*k*transpose(R_rel)*T_1
    S_el[imult_theta,itheta_2] = -alpha_stiff*k*T_2


    if mode != "driven"
        res_el[i_int] =  -k*dotp(rtmult_x,axis)
        S_el[i_int,itheta_1] = alpha_stiff*k*transpose(mult_x.v)*tilde(axis).mat*T_1
        S_el[imult_x,i_int] = alpha_stiff*k*rot(psi_1,axis).v
        S_el[i_int,imult_x] =  S_el[imult_x,i_int]'
        S_el[itheta_1,i_int] = S_el[i_int,itheta_1]'
    end
 
    if mode == "force"
        force = time_function(itime,h,params)
        p_el[i_int] = force
        res_el[i_int] += force
        force  = 0.5*(force + time_function(itime-1,h,params))
    end

    if mode == "spring"
        stiffness = params[1]
        damping = params[2]
        d_0 = params[3]
        res_el[i_int] -= (stiffness*(elong - d_0) + damping*xdot_int)
        S_el[i_int,i_int] += alpha_stiff*stiffness + theta_p*damping
    end


    if matrix == false 
        push_element_sparse(res,p,iel2,inv_loc,S_el,res_el,p_el)

        if mode == "driven"
            ext_work_el =  k*delong*transpose((R_1*axis).v)*mult_x_av
        end
        if mode == "force"
            delong = Ds[i_int]
            ext_work_el =  force*delong
        end

        mode == "spring" ? str_el = 0.5*stiffness*(elong - d_0)^2 : str_el = 0.0

        return str_el,ext_work_el
    else 
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc,S_el,p_el)
    end

end

"""
    prismatic_joint_force
"""
function prismatic_joint_force(nbr::Int,Dy::Vector,y_n::Vector,ydot::Vector,res::Vector,p::Vector,itime::Int,h::Float64)


    pjc = SetElements.prism_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, pjc.number)[1]

    str_time_function = pjc.time_function[iel]

    params = pjc.params[iel]
    k = pjc.scale_factor[iel]
    RV_rel = pjc.RV_rel[iel]
    RV_line = pjc.RV_line[iel]
    l_0 =  pjc.l_0[iel]
    
    mode = pjc.mode[iel]
    nodes =  pjc.node_orders[iel]

    axis = rot(RV_line,Vec3(1.0,0.0,0.0))
    

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    mode == "driven" ? ndim = 18 : ndim =19

    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    nodes = pjc.node_orders[iel]
    inode1 = nodes[1]
    inode2 = nodes[2]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult= ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]

    inv_loc = [inv_loc_x; inv_loc_s]

    str_time_function != " " &&  (disp = input_functions[str_time_function])

    s, Ds, sdot = pull_vectors(inv_loc, y_n, Dy, ydot)

    i_x = [i for i=1:12]
    i_x1 = i_x[1:3]
    i_x2 = i_x[7:9]
    itheta_1 = i_x[4:6]
    itheta_2 = i_x[10:12]

    if mode != "driven"
        i_int = 13
        imult= [i for i=14:19]
        x_int_n = s[i_int]
        xdot_int = sdot[i_int]
        Dx_int = Ds[i_int]
    else
        imult = [i for i=13:18]
    end

    imult_x = imult[1:3]
    imult_theta = imult[4:6]
    mult_x = Vec3(Ds[imult_x])
    mult_x_n = s[imult_x]
    mult_x_av = 0.5*(mult_x.v + mult_x_n)
    mult_theta = Vec3(Ds[imult_theta])

    Dpsi_1 = RV3(Ds[itheta_1])
    Dpsi_2 = RV3(Ds[itheta_2])
    x_1 = cf[inode1].x + Vec3(Ds[i_x1])
    x_2 = cf[inode2].x + Vec3(Ds[i_x2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)

    R_1 = (rot(psi_1)).mat

    u_psi = RV3(-RV_rel,RV3(-psi_1, psi_2))

    if mode == "driven"
        elong  = time_function(itime,h,params)
        delong  = elong - time_function(itime-1,h,params)
    else
        elong = x_int_n + Dx_int
        delong = Dx_int
    end

    u_psi = RV3(-RV_rel,RV3(-psi_1, psi_2))
    X = (l_0+elong)*axis  
    u_x = x_2 - x_1 - rot(psi_1,X)

    rtmult_x = rot(-psi_1,mult_x)

    res_el[imult_x] = k*u_x.v
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x1] = -k*mult_x.v
    res_el[i_x2] = k*mult_x.v
    res_el[itheta_1] =  - k*rot(RV_rel,mult_theta).v - k*(crossp(X,rtmult_x)).v
    res_el[itheta_2] =    k*mult_theta.v

    if mode != "driven"
        res_el[i_int] =  -k*dotp(rtmult_x,axis)
    end
 
    if mode == "force"
        force = time_function(itime,h,params)
        p_el[i_int] = force
        res_el[i_int] += force
        force  = 0.5*(force + time_function(itime-1,h,params))
    end

    if mode == "spring"
        stiffness = params[1]
        damping = params[2]
        d_0 = params[3]
        res_el[i_int] -= (stiffness*(elong - d_0) + damping*xdot_int)
    end

    for i = 1:ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    
    mode == "driven" &&  (ext_work_el = k*delong*transpose((R_1*axis).v)*mult_x_av)
    mode == "force" &&   ( ext_work_el =  force*delong)
    mode == "spring" ? str_el = 0.5*stiffness*(elong - d_0)^2 : str_el = 0.0

    
    return str_el,ext_work_el


end