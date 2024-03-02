"""
    spherical_joint
"""
function spherical_joint(nbr::Int,Dy::Vector{Float64},res::Vector{Float64},matrix::Bool)

    ndim = 15
    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))

    sjc = Main.SetElements.spherical_joint_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, sjc.numbers)[1]

    nodes = sjc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]

    relative_positions = sjc.relative_positions[iel]
    k = sjc.scale_factor[iel]

    pos1 = relative_positions[1]
    pos2 = relative_positions[2]
  

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc= [inv_loc_x; inv_loc_mult]

    s = pull_vectors(inv_loc,Dy)

    imult_x = [i for i=13:15]
    ix_1 =  [i for i=1:3]
    itheta_1 = ix_1 .+3
    ix_2 =  ix_1 .+ 6
    itheta_2 = itheta_1 .+6

    mult_x = Vec3(s[imult_x])
    Dpsi_1 = RV3(s[itheta_1])
    Dpsi_2 = RV3(s[itheta_2])
    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_2 = cf[inode2].x + Vec3(s[ix_2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)
   

    T_1 = (tang(Dpsi_1)).mat
    T_2 = (tang(Dpsi_2)).mat

    u_x = x_2 + rot(psi_2,pos2) - x_1  - rot(psi_1,pos1)

    rtmult_2 = rot(-psi_2,mult_x)
    rtmult_1 = rot(-psi_1,mult_x)

    res_el[imult_x]  = k*u_x.v
    res_el[ix_1] = - k*mult_x.v
    res_el[ix_2] =  k*mult_x.v
    res_el[itheta_2] =   + k*(crossp(pos2,rtmult_2)).v
    res_el[itheta_1] =   - k*(crossp(pos1,rtmult_1)).v 


    S_el[imult_x,ix_1] = k*eye(3)
    S_el[imult_x,ix_2] = - k*eye(3)
    S_el[ix_1,imult_x] = k*eye(3)
    S_el[ix_2,imult_x] = - k*eye(3)
    S_el[itheta_1,itheta_1] =   ((k*tilde(pos1)*tilde(rtmult_1)).mat)*T_1
    S_el[itheta_2,itheta_2] =  - ((k*tilde(pos2)*tilde(rtmult_2)).mat)*T_2
    A_1 =  (tilde(pos1)*rot(-psi_1)).mat
    A_2 =  (tilde(pos2)*rot(-psi_2)).mat
    S_el[itheta_1,imult_x] =  k*A_1
    S_el[itheta_2,imult_x] =  -k*A_2
    S_el[imult_x,itheta_1] =  k*transpose(A_1)*T_1
    S_el[imult_x,itheta_2] =  -k*transpose(A_2)*T_2


    if matrix == false

        push_element_sparse(res,iel2,inv_loc,S_el,res_el)
    else 
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end
end

"""
    spherical_joint_force
"""
function spherical_joint_force(nbr::Int,Dy::Vector{Float64},res::Vector{Float64})

    Ndim = 15
    res_el  = Vector{Float64}(zeros(Ndim))

    sjc = Main.SetElements.spherical_joint_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, sjc.numbers)[1]

    nodes = sjc.node_orders[iel]

    inode1 = nodes[1]
    inode2 = nodes[2]


    relative_positions = sjc.relative_positions[iel]
    k = sjc.scale_factor[iel]

    pos1 = relative_positions[1]
    pos2 = relative_positions[2]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc= [inv_loc_x; inv_loc_mult]

    s = pull_vectors(inv_loc,Dy)

    imult_x = [i for i=13:15]
    ix_1 =  [i for i=1:3]
    itheta_1 = ix_1 .+3
    ix_2 =  ix_1 .+ 6
    itheta_2 = itheta_1 .+6

    mult_x = Vec3(s[imult_x])
    Dpsi_1 = RV3(s[itheta_1])
    Dpsi_2 = RV3(s[itheta_2])
    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_2 = cf[inode2].x + Vec3(s[ix_2])
    psi_1 = RV3(cf[inode1].psi,Dpsi_1)
    psi_2 = RV3(cf[inode2].psi,Dpsi_2)
   
    u_x = x_2 + rot(psi_2,pos2) - x_1  - rot(psi_1,pos1)

    u_x = x_2 + rot(psi_2,pos2) - x_1  - rot(psi_1,pos1)

    rtmult_2 = rot(-psi_2,mult_x)
    rtmult_1 = rot(-psi_1,mult_x)

    res_el[imult_x]  = k*u_x.v
    res_el[ix_1] = - k*mult_x.v
    res_el[ix_2] =  k*mult_x.v
    res_el[itheta_2] =   + k*(crossp(pos2,rtmult_2)).v
    res_el[itheta_1] =   - k*(crossp(pos1,rtmult_1)).v 

    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i])
    end

    return

end

