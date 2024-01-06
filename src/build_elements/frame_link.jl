"""
    frame_link
"""
function frame_link(nbr::Int,Dy::Vector{Float64},res::Vector{Float64},matrix::Bool)

    S_el = Array{Float64,2}(zeros(18,18))
    res_el  = Vector{Float64}(zeros(18))

    flc = Main.SetElements.frame_link_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, flc.numbers)[1]

    nodes = flc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]


    X_0 = flc.relative_position[iel]
    RV_0 = flc.relative_orientation[iel]
    k = flc.scale_factor[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc= [inv_loc_x; inv_loc_mult]

    s = pull_vectors(inv_loc,Dy)

    imult_x = [i for i=13:15]
    imult_theta = imult_x .+ 3
    ix_1 =  [i for i=1:3]
    itheta_1 = ix_1 .+3
    ix_2 =  ix_1 .+ 6
    itheta_2 = itheta_1 .+6

    mult_x = Vec3(s[imult_x])
    mult_theta = s[imult_theta]
    Dpsi_1 = RV3(s[itheta_1])
    Dpsi_2 = RV3(s[itheta_2])
    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_2 = cf[inode2].x + Vec3(s[ix_2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)
   

    T_1 = (tang(Dpsi_1)).mat
    T_2 = (tang(Dpsi_2)).mat
    R_0 = (rot(RV_0)).mat

    u_psi = RV3(-RV_0,RV3(-psi_1, psi_2))
    u_x = x_2 - x_1  - rot(psi_1,X_0)

    rtmult_x=rot(-psi_1,mult_x)

    res_el[imult_x]  = k*u_x.v
    res_el[imult_theta] = u_psi.v
    res_el[ix_1] = - k*mult_x.v
    res_el[ix_2] =  k*mult_x.v
    res_el[itheta_1] =  - k*R_0*mult_theta - k*(crossp(X_0,rtmult_x)).v
    res_el[itheta_2] =    k*mult_theta


    S_el[imult_x,ix_1] = k*eye(3)
    S_el[imult_x,ix_2] = - k*eye(3)
    S_el[ix_1,imult_x] = k*eye(3)
    S_el[ix_2,imult_x] = - k*eye(3)
    S_el[itheta_1,itheta_1] =   ((k*tilde(X_0)*tilde(rtmult_x)).mat)*T_1
    A =  (tilde(X_0)*rot(-psi_1)).mat
    S_el[itheta_1,imult_x] =  k*A
    S_el[imult_x,itheta_1] =  k*transpose(A)*T_1

    S_el[imult_theta,itheta_1] = k*transpose(R_0)*T_1
    S_el[imult_theta,itheta_2] = -k*T_2
    S_el[itheta_1,imult_theta] = k*R_0
    S_el[itheta_2,imult_theta] = -k*eye(3)


    if matrix == false

        push_element_sparse(res,iel2,inv_loc,S_el,res_el)
    else 
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end
end

"""
    frame_link_force
"""
function frame_link_force(nbr::Int,Dy::Vector{Float64},res::Vector{Float64})

    Ndim = 18
    res_el  = Vector{Float64}(zeros(Ndim))

    flc = Main.SetElements.frame_link_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, flc.numbers)[1]

    nodes = flc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]


    X_0 = flc.relative_position[iel]
    RV_0 = flc.relative_orientation[iel]
    k = flc.scale_factor[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc= [inv_loc_x; inv_loc_mult]

    s = pull_vectors(inv_loc,Dy)

    imult_x = [i for i=13:15]
    imult_theta = imult_x .+ 3
    ix_1 =  [i for i=1:3]
    itheta_1 = ix_1 .+3
    ix_2 =  ix_1 .+ 6
    itheta_2 = itheta_1 .+6

    mult_x = Vec3(s[imult_x])
    mult_theta = s[imult_theta]
    Dpsi_1 = RV3(s[itheta_1])
    Dpsi_2 = RV3(s[itheta_2])
    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_2 = cf[inode2].x + Vec3(s[ix_2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)
   

    R_0 = (rot(RV_0)).mat

    u_psi = RV3(-RV_0,RV3(-psi_1, psi_2))
    u_x = x_2 - x_1  - rot(psi_1,X_0)

    rtmult_x=rot(-psi_1,mult_x)

    res_el[imult_x]  = k*u_x.v
    res_el[imult_theta] = u_psi.v
    res_el[ix_1] = - k*mult_x.v
    res_el[ix_2] =  k*mult_x.v
    res_el[itheta_1] =  - k*R_0*mult_theta - k*(crossp(X_0,rtmult_x)).v
    res_el[itheta_2] =    k*mult_theta

    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i])
    end

    return

end

