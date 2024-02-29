"""
    hinge
"""
function hinge(nbr::Int,Dy::Vector,y_n::Vector,res::Vector,p::Vector,niter::Int,itime::Int,h::Float64,matrix::Bool)


    hc = Main.SetElements.hinge_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions

    iel = findfirst(x -> x == nbr, hc.numbers)[1]
    mode = hc.mode[iel]
    str_time_function = hc.time_function[iel]
    k = hc.scale_factor[iel]
    params = hc.params[iel]

    time_function = input_functions[str_time_function]

    mode == "driven" ? ndim = 18 : ndim =19


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    nodes = hc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]

    X1 = hc.positions[iel][1]
    X2 = hc.positions[iel][2]
    or1 = hc.orientations[iel][1]
    or2 = hc.orientations[iel][2]
    axis = hc.axis[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]

    inv_loc = [inv_loc_x; inv_loc_s]

    #
    # The hinge angle is a variable with no dynamics.
    # Its value reults from the expression of the constraints.
    # It must thus be restset to zero at the beginning of each time step.
    #


    if  mode != "driven"
        s = pull_vectors(inv_loc, Dy)
        if niter == 1
             Dy[inv_loc_int] .= 0.0
             y_n[inv_loc_int] .= 0.0
        end
    else
        s_n, s = pull_vectors(inv_loc, y_n, Dy)
    end

    i_x = [i for i=1:12]
    i_x1 = i_x[1:3]
    i_x2 = i_x[7:9]
    i_theta1 = i_x[4:6]
    i_theta2 = i_x[10:12]

    if mode != "driven"
        i_int = 13
        imult= [i for i=14:19]
        x_int = s[i_int]
        mult = s[imult]
    else
        imult = [i for i=13:18]
        mult = s[imult]
        mult_theta_av = 0.5*(s_n[imult[4:6]]+s[imult[4:6]])
    end

    imult_x = imult[1:3]
    imult_theta = imult[4:6]
    mult_x = Vec3(mult[1:3])
    mult_theta = Vec3(mult[4:6])

    Dpsi_1 = RV3(s[i_theta1])
    Dpsi_2 = RV3(s[i_theta2])
    x1 = cf[inode1].x + Vec3(s[i_x1])
    x2 = cf[inode2].x + Vec3(s[i_x2])
    psi1 = RV3(cf[inode1].psi,Dpsi_1)
    psi2 = RV3(cf[inode2].psi,Dpsi_2)


    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        psi_hinge = RV3(x_int*xi)
    end



    R1 = rot(psi1)
    R2 = rot(psi2)
    T1 = tang(Dpsi_1)
    T2 = tang(Dpsi_2)

    R1r = rot(or1)
    R2r = rot(or2)

    u_x = x2 + rot(psi2,X2) - x1 - rot(psi1,X1)

    psi_rel = RV3(-RV3(psi1,or1), RV3(psi2,or2))
    R_hinge = rot(psi_hinge)
    u_psi = RV3(-psi_hinge, psi_rel)

    res_el[imult_x] = k*u_x.v
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x1] = - k*mult_x.v
    res_el[i_x2] =   k*mult_x.v
    res_el[i_theta1] =  -k*(crossp(X1,rot(-psi1,mult_x))  +  rot(RV3(or1,psi_hinge),mult_theta)).v
    res_el[i_theta2] =   k*(crossp(X2,rot(-psi2,mult_x)) + rot(or2,mult_theta)).v

    if mode != "driven"
        res_el[i_int] =  -k*transpose(xi)*mult_theta.v
    end

    S_el[imult_x,i_x1] = k*eye(3)
    S_el[imult_x,i_x2] = - k*eye(3)
    S_el[imult_x,i_theta1] = -k*(R1*tilde(X1)*T1).mat
    S_el[imult_x,i_theta2] =  k*(R2*tilde(X2)*T2).mat
    S_el[i_x1,imult_x] = k*eye(3)
    S_el[i_x2,imult_x] = -k*eye(3)
    S_el[i_theta1,imult_x] = -k*transpose(R1*tilde(X1)).mat
    S_el[i_theta2,imult_x] =  k*transpose(R2*tilde(X2)).mat
    S_el[i_theta1,i_theta1]  =  k*(tilde(X1)*tilde(transpose(R1)*mult_x)*T1).mat
    S_el[i_theta2,i_theta2]  = - k*(tilde(X2)*tilde(transpose(R2)*mult_x)*T2).mat

    S_el[imult_theta,i_theta1] =    k*(transpose(R1r*R_hinge)*T1).mat
    S_el[imult_theta,i_theta2] =  - k*(transpose(R2r)*T2).mat
    S_el[i_theta1,imult_theta] =    k*(R1r*R_hinge).mat
    S_el[i_theta2,imult_theta] =  - k*R2r.mat


    if mode != "driven"
        S_el[imult_theta,i_int] = k*xi
        S_el[i_int,imult_theta] = k*transpose(xi)
    end

    if mode == "force"
        torque = time_function(itime,h,params)
        p_el[i_int] = torque
        res_el[i_int] += torque
        torque = 0.5*(torque+time_function(itime-1,h,params))
#        S_el[i_theta1, i_int] = k*(rot(RV3(or1,psi_hinge))*tilde(axis)).mat*mult_theta
    end

    if matrix == false

        push_element_sparse(res,p,iel2,inv_loc,S_el,res_el,p_el)

        if mode == "driven"
            ext_work_el =  k*alpha*transpose(xi)*mult_theta_av
        end
        if mode == "force"
            phi = Dy[inv_loc_int][1]
            ext_work_el =  torque*phi
        end

        return ext_work_el
    else 
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end

end

"""
    hinge_force
"""
function hinge_force(nbr::Int,Dy::Vector,y_n::Vector,res::Vector,p::Vector,niter::Int,itime::Int,h::Float64)


    hc = Main.SetElements.hinge_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions
    
    iel = findfirst(x -> x == nbr, hc.numbers)[1]
    mode = hc.mode[iel]
    str_time_function = hc.time_function[iel]
    k = hc.scale_factor[iel]
    params = hc.params[iel]

    time_function = input_functions[str_time_function]

    mode == "driven" ? ndim = 18 : ndim =19


    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    nodes = hc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]

    X1 = hc.positions[iel][1]
    X2 = hc.positions[iel][2]
    or1 = hc.orientations[iel][1]
    or2 = hc.orientations[iel][2]
    axis = hc.axis[iel]


    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]

    inv_loc = [inv_loc_x; inv_loc_s]

    #
    # The hinge angle is a variable with no dynamics.
    # Its value reults from the expression of the constraints.
    # It must thus be restset to zero at the beginning of each time step.
    #


    if  mode != "driven"
        s = pull_vectors(inv_loc, Dy)
        if niter == 1
             Dy[inv_loc_int] .= 0.0
             y_n[inv_loc_int] .= 0.0
        end
    else
        s_n, s = pull_vectors(inv_loc, y_n, Dy)
    end

    i_x = [i for i=1:12]
    i_x1 = i_x[1:3]
    i_x2 = i_x[7:9]
    i_theta1 = i_x[4:6]
    i_theta2 = i_x[10:12]

    if mode != "driven"
        i_int = 13
        imult= [i for i=14:19]
        x_int = s[i_int]
        mult = s[imult]
    else
        imult = [i for i=13:18]
        mult = s[imult]
        mult_theta_av = 0.5*(s_n[imult[4:6]]+s[imult[4:6]])
    end

    imult_x = imult[1:3]
    imult_theta = imult[4:6]
    mult_x = Vec3(mult[1:3])
    mult_theta = Vec3(mult[4:6])


    Dpsi_1 = RV3(s[i_theta1])
    Dpsi_2 = RV3(s[i_theta2])
    x1 = cf[inode1].x + Vec3(s[i_x1])
    x2 = cf[inode2].x + Vec3(s[i_x2])
    psi1 = RV3(cf[inode1].psi,Dpsi_1)
    psi2 = RV3(cf[inode2].psi,Dpsi_2)


    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        psi_hinge = RV3(x_int*xi)
    end

    u_x = x2 + rot(psi2,X2) - x1 - rot(psi1,X1)

    psi_rel = RV3(-RV3(psi1,or1), RV3(psi2,or2))
    u_psi = RV3(-psi_hinge, psi_rel)

    res_el[imult_x] = k*u_x.v
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x1] = - k*mult_x.v
    res_el[i_x2] =   k*mult_x.v
    res_el[i_theta1] =  -k*(crossp(X1,rot(-psi1,mult_x))  +  rot(RV3(or1,psi_hinge),mult_theta)).v
    res_el[i_theta2] =   k*(crossp(X2,rot(-psi2,mult_x)) + rot(or2,mult_theta)).v

    if mode != "driven"
        res_el[i_int] =  -k*transpose(xi)*mult_theta.v
    end

    if mode == "force"
        torque = time_function(itime,h,params)
        p_el[i_int] = torque
        res_el[i_int] += torque
        torque = 0.5*(torque+time_function(itime-1,h,params))

    end

    for i = 1:ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    if mode == "driven"
        ext_work_el =  k*alpha*transpose(xi)*mult_theta_av
    end
    if mode == "force"
        phi = Dy[inv_loc_int][1]
        ext_work_el =  torque*phi
    end

    return ext_work_el

end

