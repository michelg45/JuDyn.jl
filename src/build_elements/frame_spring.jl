"""
    frame_spring
"""
function frame_spring(nbr::Int,res::Vector)

    S_el = Array{Float64,2}(zeros(12,12))
    res_el  = Vector{Float64}(zeros(12))

    fsc = SetElements.frame_spring_container
    nbc = node_container
    cf =  Frames.current_frames
    ec =  element_container


    iel = findfirst(x -> x == nbr, fsc.numbers)[1]

    nodes = fsc.node_orders[iel]

    inode1 = findfirst(x -> x == nodes[1], nbc.node_numbers)[1]
    inode2 = findfirst(x -> x == nodes[2], nbc.node_numbers)[1]

    k_x = fsc.extension_stiffness[iel]
    k_theta = fsc.rotation_stiffness[iel]
    X_0 = fsc.initial_displacements[iel]
    RV_0 = fsc.initial_rotations[iel]
    K = diagonal([k_x.v;k_theta.v])

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc = ec.inv_loc_x[iel2]

    x  =  Vector{Vec3}([cf[inode1].x, cf[inode2].x])
    psi  = Vector{RV3}([cf[inode1].psi, cf[inode2].psi])
    psi[2] = RV3(-RV_0,psi[2])
    u_psi = RV3(-psi[1], psi[2])
    PSI = norm(u_psi.v)
    u_x = invtang(-u_psi)*(rot(-psi[1],x[2]- x[1])-X_0)


    epsilon = [u_x.v ; u_psi.v]

    minvtang_1 = -invtang(-u_psi).mat
    invtang_2 = invtang(u_psi).mat
    mDinvtang_1 = - Dinvtang(-u_x,-u_psi).mat
    Dinvtang_2 = Dinvtang(u_x,u_psi).mat
    meye3 = -eye(3)
    P_L = [ minvtang_1  mDinvtang_1 invtang_2  Dinvtang_2; zeros(3,3) meye3 zeros(3,3) eye(3)]
    P_R = [ minvtang_1  mDinvtang_1 invtang_2  Dinvtang_2; zeros(3,3)  minvtang_1 zeros(3,3)  invtang_2]
    Keps = K*epsilon
    res_el = - transpose(P_L)*Keps
    S_el = transpose(P_L)*K*P_R
    str_el = 0.5*transpose(epsilon)*Keps

    push_element_sparse(res,iel2,inv_loc,S_el,res_el)

    return str_el


end
