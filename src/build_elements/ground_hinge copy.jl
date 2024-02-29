"""
    ground_hinge
"""
function ground_hinge(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},res::Vector{Float64},p::Vector{Float64},niter::Int,itime::Int,h::Float64,matrix::Bool)


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

    time_function = input_functions[str_time_function]

    mode == "driven" ? ndim = 6 : ndim =7


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    inode = ghc.node_orders[iel]

    # inode = findfirst(x -> x == node, nbc.node_numbers)[1]

    RV_0 = ghc.orientation[iel]
    axis = ghc.axis[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_theta= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]
    inv_loc = [inv_loc_theta; inv_loc_s]

    #
    # The hinge angle is a variable with no dynamics.
    # Its value reults from the expression of the constraints.
    # It must thus be restset to zero at the beginning of each time step.
    #

    if  mode != "driven" && niter == 1
         Dy[inv_loc_int] .= 0.0
    end    

    s, Ds = pull_vectors(inv_loc, y_n, Dy)




    i_theta = [i for i=1:3]

    if mode != "driven"
        i_int = 4
        imult_theta = i_theta .+ 4
        x_int_n = s[i_int]
        Dx_int = Ds[i_int]
        mult_theta = Ds[imult_theta]
    else
        imult_theta = i_theta .+ 3
        mult_theta = Ds[imult_theta]
        mult_theta_av = 0.5*(s[imult_theta]+Ds[imult_theta])
    end


    
    Dpsi  = RV3(Ds[i_theta])
    psi   = RV3(cf[inode].psi,Dpsi)


    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        #x_int = x_int_n + Dx_int
        x_int = Dx_int
        psi_hinge = RV3(x_int*xi)
    end

    T =  tang(Dpsi).mat
    
    u_psi = RV3(-RV_0,RV3(-psi_hinge, psi))

    res_el[imult_theta] = k*u_psi.v
    res_el[i_theta] =    k*mult_theta
    S_el[imult_theta,i_theta] = -k*T
    S_el[i_theta,imult_theta] = - k*eye(3)

    if mode != "driven"
        res_el[i_int] =    -k*transpose(xi)*mult_theta
        S_el[imult_theta,i_int] = k*xi
        S_el[i_int,imult_theta] = k*transpose(xi)
    end

    if mode == "force"
        torque = time_function(itime,h,params)
        p_el[i_int] = torque
        res_el[i_int] += torque
        torque = 0.5*(torque+time_function(itime-1,h,params))
    end

    if matrix == false
        push_element_sparse(res,p,iel2,inv_loc,S_el,res_el,p_el)

        mode == "driven" &&  (ext_work_el =  k*alpha*transpose(xi)*mult_theta_av)
        mode == "force" &&   ( ext_work_el =  (torque*Dx_int)[1])
        return ext_work_el
    else
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end

end

"""
    ground_hinge_force
"""
function ground_hinge_force(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},res::Vector{Float64},p::Vector{Float64},niter::Int,itime::Int,h::Float64)


    ghc = SetElements.ground_hinge_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, ghc.numbers)[1]
    mode = ghc.mode[iel]
    str_time_function = ghc.time_function[iel]
    params = ghc.params[iel]
    k = ghc.scale_factor[iel]

    time_function = input_functions[str_time_function]

    mode == "driven" ? ndim = 6 : ndim =7


    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    inode = ghc.node_orders[iel]

    # inode = findfirst(x -> x == node, nbc.node_numbers)[1]

    RV_0 = ghc.orientation[iel]
    axis = ghc.axis[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_theta= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]
    inv_loc = [inv_loc_theta; inv_loc_s]

    #
    # The hinge angle is a variable with no dynamics.
    # Its value reults from the expression of the constraints.
    # It must thus be restset to zero at the beginning of each time step.
    #
    if  mode != "driven" && niter == 1
         Dy[inv_loc_int] .= 0.0
    end

    s, Ds = pull_vectors(inv_loc, y_n, Dy)

    i_theta = [i for i=1:3]

    if mode != "driven"
        i_int = 4
        imult_theta = i_theta .+ 4
        x_int_n = s[i_int]
        Dx_int = Ds[i_int]
        mult_theta = Ds[imult_theta]
    else
        imult_theta = i_theta .+ 3
        mult_theta = Ds[imult_theta]
        mult_theta_av = 0.5*(s[imult_theta]+Ds[imult_theta])
    end

    Dpsi  = RV3(Ds[i_theta])
    psi   = RV3(cf[inode].psi,Dpsi)

    xi = axis.v

    if mode == "driven"
        alpha = time_function(itime,h,params)
        psi_hinge = RV3(alpha*xi)
        alpha -= time_function(itime-1,h,params)
    else
        x_int = Dx_int
        psi_hinge = RV3(x_int*xi)
    end

    u_psi = RV3(-RV_0,RV3(-psi_hinge, psi))



    res_el[imult_theta] = k*u_psi.v
    res_el[i_theta] =    k*mult_theta

    if mode != "driven"
        res_el[i_int] =    -k*transpose(xi)*mult_theta
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

    mode == "driven" &&  (ext_work_el =  k*alpha*transpose(xi)*mult_theta_av)
    mode == "force" &&   ( ext_work_el =  (torque*Dx_int)[1])

    return ext_work_el

end

