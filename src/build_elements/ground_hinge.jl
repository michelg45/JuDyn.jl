"""
    ground_hinge
"""
function ground_hinge(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},ydot::Vector{Float64},res::Vector{Float64},p::Vector{Float64},alpha_stiff::Float64,theta_p::Float64,matrix::Bool,niter::Int,itime::Int,h::Float64)


    ghc = SetElements.ground_hinge_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, ghc.numbers)[1]
    mode = ghc.mode[iel]
    str_time_function = ghc.time_function[iel]
    params = ghc.params[iel]
    k = ghc.scale_factor[iel]

    str_time_function != " " &&  (time_function = input_functions[str_time_function])

    mode == "driven" ? ndim = 12 : ndim =13


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    inode = ghc.node_orders[iel]

    # inode = findfirst(x -> x == node, nbc.node_numbers)[1]

    RV_0 = ghc.orientation[iel]
    X_0 = ghc.position[iel]
    axis = ghc.axis[iel]

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

    s, Ds, sdot = pull_vectors(inv_loc, y_n, Dy, ydot)

    i_x = [i for i=1:3]
    i_theta = i_x .+ 3

    if mode != "driven"
        i_int = 7
        imult_x = i_x .+ 7
        imult_theta = i_theta .+ 7
        x_int = s[i_int]
        xdot_int = sdot[i_int]
        Dx_int = Ds[i_int]
        x_int += Dx_int
        mult_x = Ds[imult_x]
        mult_theta = Ds[imult_theta]
    else
        imult_x = i_x .+ 6
        imult_theta = i_theta .+ 6
        mult_x = Ds[imult_x]
        mult_theta = Ds[imult_theta]
        mult_theta_av = 0.5*(s[imult_theta]+Ds[imult_theta])
    end

    x = (cf[inode].x).v  + Ds[i_x]
    Dpsi  = RV3(Ds[i_theta])
    psi   = RV3(cf[inode].psi,Dpsi)


    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        psi_hinge = RV3(x_int*xi)
    end

    

    T =  tang(Dpsi).mat
    
    u_x = x - X_0.v
    u_psi = RV3(-psi_hinge,RV3(-RV_0,psi))
    res_el[imult_x] = k*u_x
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x] = k*mult_x
    res_el[i_theta] = k*mult_theta
    S_el[imult_x,i_x] = -alpha_stiff*k*eye(3)
    S_el[i_x,imult_x] = -alpha_stiff*k*eye(3)
    S_el[imult_theta,i_theta] = -alpha_stiff*k*T
    S_el[i_theta,imult_theta] = - alpha_stiff*k*eye(3)

    if mode != "driven"
        res_el[i_int] = -k*transpose(xi)*mult_theta
        S_el[i_int,imult_theta] = k*transpose(xi)
        S_el[imult_theta,i_int] = alpha_stiff*k*xi
        S_el[i_int,imult_theta] = alpha_stiff*k*transpose(xi)
    end

    if mode == "force"
        torque = time_function(itime,h,params)
        p_el[i_int] = torque
        res_el[i_int] += torque
        torque = 0.5*(torque+time_function(itime-1,h,params))
    end

    if mode == "spring"
        stiffness = params[1]
        damping = params[2]
        alpha_0 = params[3]
        res_el[i_int] -= (stiffness*(x_int - alpha_0) + damping*xdot_int)
        S_el[i_int,i_int] += alpha_stiff*stiffness + theta_p*damping
    end

    mode == "spring" ? str_el = 0.5*stiffness*(x_int - alpha_0)^2 : str_el = 0.0

    if matrix == false
        push_element_sparse(res,p,iel2,inv_loc,S_el,res_el,p_el)

        mode == "driven" &&  (ext_work_el =  k*alpha*transpose(xi)*mult_theta_av)
        mode == "force" &&   ( ext_work_el =  (torque*Dx_int)[1])
        return str_el, ext_work_el
    else
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end

end

"""
    ground_hinge_force
"""
function ground_hinge_force(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},ydot::Vector{Float64},res::Vector{Float64},p::Vector{Float64},niter::Int,itime::Int,h::Float64)


    ghc = SetElements.ground_hinge_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, ghc.numbers)[1]
    mode = ghc.mode[iel]
    str_time_function = ghc.time_function[iel]
    params = ghc.params[iel]
    k = ghc.scale_factor[iel]

    str_time_function != " " &&  (time_function = input_functions[str_time_function])

    mode == "driven" ? ndim = 12 : ndim =13

    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    inode = ghc.node_orders[iel]

    # inode = findfirst(x -> x == node, nbc.node_numbers)[1]

    RV_0 = ghc.orientation[iel]
    X_0 = ghc.position[iel]
    axis = ghc.axis[iel]

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

    s, Ds, sdot = pull_vectors(inv_loc, y_n, Dy, ydot)

    i_x = [i for i=1:3]
    i_theta = i_x .+ 3

    if mode != "driven"
        i_int = 7
        imult_x = i_x .+ 7
        imult_theta = i_theta .+ 7
        x_int = s[i_int]
        xdot_int = sdot[i_int]
        Dx_int = Ds[i_int]
        x_int += Dx_int
        mult_x = Ds[imult_x]
        mult_theta = Ds[imult_theta]
    else
        imult_x = i_x .+ 6
        imult_theta = i_theta .+ 6
        mult_x = Ds[imult_x]
        mult_theta = Ds[imult_theta]
        mult_theta_av = 0.5*(s[imult_theta]+Ds[imult_theta])
    end

    x = (cf[inode].x).v  + Ds[i_x]
    Dpsi  = RV3(Ds[i_theta])
    psi   = RV3(cf[inode].psi,Dpsi)


    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        psi_hinge = RV3(x_int*xi)
    end

    

    T =  tang(Dpsi).mat
    
    u_x = x - X_0.v
    u_psi = RV3(-psi_hinge,RV3(-RV_0,psi))
    res_el[imult_x] = k*u_x
    res_el[imult_theta] = k*u_psi.v
    res_el[i_x] = k*mult_x
    res_el[i_theta] = k*mult_theta

    if mode != "driven"
        res_el[i_int] = -k*transpose(xi)*mult_theta
    end

    if mode == "force"
        torque = time_function(itime,h,params)
        p_el[i_int] = torque
        res_el[i_int] += torque
        torque = 0.5*(torque+time_function(itime-1,h,params))
    end

    if mode == "spring"
        stiffness = params[1]
        damping = params[2]
        alpha_0 = params[3]
        res_el[i_int] -= (stiffness*(x_int - alpha_0) + damping*xdot_int)
    end

    for i = 1:ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    mode == "driven" &&  (ext_work_el =  k*alpha*transpose(xi)*mult_theta_av)
    mode == "force" &&   ( ext_work_el =  (torque*Dx_int)[1])
    mode == "spring" ? str_el = 0.5*stiffness*(x_int - alpha_0)^2 : str_el = 0.0

    return str_el, ext_work_el

end

